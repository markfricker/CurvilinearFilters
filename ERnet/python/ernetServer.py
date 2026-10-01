#!/usr/bin/env python
"""ernetServer.py  —  persistent ERnet worker for MATLAB integration.

Start once from PowerShell or via MATLAB (Task Scheduler); handles all ERnet
inference requests for the session without per-call model-load overhead.
Same request/response protocol as nerdyServer.py.

Usage:
    python ernetServer.py <watch_dir> [model_dir]

<watch_dir>  — directory watched for *.req.mat files written by ernetEnhance.m
<model_dir>  — folder holding the two files setupPythonEnvs('ernet') downloads
               from the ERnet-v2 GitHub repository (GPL-3.0, not bundled).
               Default: <venv root>/ernet, i.e. %USERPROFILE%\venvs\nerdy\ernet
               when run by the nERdy venv's python.exe (kept implicit so the
               Task Scheduler /TR command stays under its 261-character limit):
                 swinir_rcab_arch.py                  (model architecture)
                 20220306_ER_4class_swinir_nch1.pth   (trained weights, v2.0 release)

Request fields (*.req.mat):
    I         — 2-D float32 image in [0, 1]
    device    — 'auto', 'cuda', or 'cpu'  (default 'cpu')
    threshold — NaN = most likely class (argmax); a value in [0, 1] is a
                fixed threshold on P(ER) = 1 - P(background)
    timeout   — seconds  (default 300)

Result file (*.res.mat):
    R — binary single [0, 1], same spatial size as I; 1 = any ER class
        (tubule, sheet or sheet-based tubule)
    L — uint8 class map, same size: 0 background, 1 tubule, 2 sheet,
        3 sheet-based tubule (SBT) -- ERnet's own argmax labels after its
        class-order fix, independent of threshold

The network is the single-frame SwinIR model (4 classes) — the only ERnet
model with released weights. Preprocessing copies ERnet-v2's
Inference/model_evaluation.py (min-max, then 1-99 percentile stretch, then
8-bit quantisation), so results match the published tool. The tubule/sheet
split depends strongly on the image scale and is not returned; the
ER/background split does not (checked 0.5x-1.5x on Airyscan ER, 2026-09-30).

Stop by creating 'exit.req' in watch_dir, or kill the process.
"""

import os, sys

# ---- Strip MATLAB runtime DLLs from PATH before any other imports -------------
# Task Scheduler / schtasks inherits MATLAB's environment.  MATLAB ships its
# own MKL and C++ runtime; Python's scipy/torch find them first and deadlock.
os.environ['PATH'] = os.pathsep.join(
    p for p in os.environ.get('PATH', '').split(os.pathsep)
    if not any(s in p.upper() for s in ('MATLAB', 'MWE_', 'POLYSPACE')))
for _v in ('PYTHONPATH', 'PYTHONHOME', 'PYTHONSTARTUP'):
    os.environ.pop(_v, None)

import time, traceback
from pathlib import Path

if len(sys.argv) < 2:
    print('Usage: ernetServer.py <watch_dir> [model_dir]', file=sys.stderr)
    sys.exit(1)

watch_dir = Path(sys.argv[1])
watch_dir.mkdir(parents=True, exist_ok=True)
if len(sys.argv) >= 3:
    model_dir = Path(sys.argv[2])
else:
    model_dir = Path(sys.executable).resolve().parent.parent / 'ernet'

WEIGHTS = '20220306_ER_4class_swinir_nch1.pth'
ARCH    = 'swinir_rcab_arch.py'
WINDOW  = 4   # SwinIR window size: H and W must be multiples of this

# ---- Write PID file immediately so MATLAB knows the process started ----------
pid_file = watch_dir / 'server.pid'
pid_file.write_text(str(os.getpid()))
print(f'[server] PID {os.getpid()}  watching {watch_dir}', flush=True)

# ---- Startup diagnostic log (check <watch_dir>/server_startup.log on failure)
_log_path = watch_dir / 'server_startup.log'
def _log(msg):
    with open(_log_path, 'a') as _f:
        _f.write(f'{time.time():.3f}  {msg}\n')
        _f.flush()

_log(f'started  PID={os.getpid()}  Python={sys.version.split()[0]}')
_log(f'model_dir={model_dir}')
_log(f'PATH={os.environ.get("PATH","")[:300]}')

_log('importing numpy...')
import numpy as np
_log('importing scipy.io...')
import scipy.io as sio
_log('basic imports done')

# ---- Model singleton -----------------------------------------------------------
_model      = None
_model_dev  = None   # the device_str string that was used to load the current model
_torch_dev  = None   # the resolved torch.device

def _get_model(device_str):
    """Load the ERnet model on first call, or when device changes."""
    global _model, _model_dev, _torch_dev

    if _model is not None and device_str == _model_dev:
        return _model, _torch_dev

    _log(f'importing torch (device={device_str})...')
    import torch
    _log(f'torch ok: {torch.__version__}')

    if device_str == 'auto':
        dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        dev = torch.device(device_str)
    _log(f'resolved device: {dev}')

    for f in (ARCH, WEIGHTS):
        if not (model_dir / f).exists():
            raise FileNotFoundError(
                f'ERnet file not found: {model_dir / f}\n'
                "Run setupPythonEnvs('ernet') in MATLAB to download it.")

    sys.path.insert(0, str(model_dir))
    _log('importing swinir_rcab_arch.SwinIR_RCAB...')
    from argparse import Namespace
    from swinir_rcab_arch import SwinIR_RCAB

    _log(f'loading weights from {model_dir / WEIGHTS}...')
    t0 = time.time()
    # Constructor arguments as in ERnet-v2 Inference/model_evaluation.py
    m  = SwinIR_RCAB(Namespace(task='segment'), img_size=128, in_chans=1,
                     upscale=1, use_checkpoint=False, vis=False)
    ck = torch.load(str(model_dir / WEIGHTS), map_location=dev, weights_only=True)
    m.load_state_dict(ck['state_dict'], strict=True)
    m  = m.to(dev)
    m.eval()
    _log(f'model ready in {time.time()-t0:.1f}s on {dev}')
    print(f'[server] model ready in {time.time()-t0:.1f}s on {dev}', flush=True)

    _model     = m
    _model_dev = device_str
    _torch_dev = dev
    return _model, _torch_dev


def _preprocess(I):
    """ERnet-v2 model_evaluation.py: min-max, 1-99 percentile stretch, 8-bit."""
    I = I.astype(np.float64)
    lo, hi = I.min(), I.max()
    if hi <= lo:
        return np.zeros(I.shape, dtype=np.float32)
    I = (I - lo) / (hi - lo)
    p1, p99 = np.percentile(I, 1), np.percentile(I, 99)
    if p99 > p1:
        I = np.clip((I - p1) / (p99 - p1), 0, 1)
    # uint8 quantisation, as ERnet's PIL / ToTensor round trip
    return (I * 255).astype(np.uint8).astype(np.float32) / 255


# ---- Pre-load on CPU so first request is fast --------------------------------
try:
    _get_model('cpu')
    (watch_dir / 'server.ready').write_text('1')
    print('[server] ready', flush=True)
    _log('server.ready written')
except Exception as _e:
    _log(f'WARNING: model pre-load failed: {_e}\n{traceback.format_exc()}')
    print(f'[server] WARNING: model pre-load failed: {_e}', flush=True)
    # Still write ready so MATLAB doesn't time out; first request will try again
    (watch_dir / 'server.ready').write_text('1')


# ---- Request loop -------------------------------------------------------------
while True:
    # Graceful exit
    if (watch_dir / 'exit.req').exists():
        try:
            (watch_dir / 'exit.req').unlink()
        except Exception:
            pass
        print('[server] exit requested — shutting down', flush=True)
        break

    for req_file in sorted(watch_dir.glob('*.req.mat')):
        stem     = req_file.stem          # e.g. '20260601_123456_42.req'
        base     = stem.replace('.req', '')
        res_file = watch_dir / f'{base}.res.mat'
        err_file = watch_dir / f'{base}.err'

        try:
            # Retry loadmat briefly — belt-and-braces guard for mid-write reads
            for _attempt in range(5):
                try:
                    mat = sio.loadmat(str(req_file))
                    break
                except Exception:
                    if _attempt == 4:
                        raise
                    time.sleep(0.1)

            def _scalar(key, default):
                if key not in mat:
                    return default
                v = np.squeeze(mat[key])
                return float(v.item()) if v.ndim == 0 else default

            I          = np.squeeze(mat['I']).astype(np.float32)
            device_str = str(mat['device'].flat[0]) if 'device' in mat else 'cpu'
            threshold  = _scalar('threshold', float('nan'))
            timeout    = _scalar('timeout', 300.0)

            m, dev = _get_model(device_str)

            _result = [None, None]

            def _run():
                try:
                    import torch

                    H, W = I.shape
                    x = _preprocess(I)
                    # Pad to a multiple of the window size (window_partition
                    # needs whole windows); 'symmetric' works for any size
                    ph, pw = (-H) % WINDOW, (-W) % WINDOW
                    x = np.pad(x, ((0, ph), (0, pw)), mode='symmetric')
                    t = torch.from_numpy(x)[None, None].to(dev)

                    with torch.no_grad():
                        out = m(t)[0]            # 4 x H' x W' class scores
                    out = out[:, :H, :W]

                    # Channel 0 is background; 1-3 are the ER classes
                    # (ERnet's own post-hoc relabelling only permutes 1-3)
                    raw = out.argmax(0)
                    if np.isnan(threshold):
                        R = (raw != 0)
                    else:
                        pEr = 1 - torch.softmax(out, 0)[0]
                        R = (pEr >= threshold)
                    # Class map in ERnet's published order. Raw channel
                    # indices map as model_evaluation.py's workaround:
                    # raw 1 -> SBT, raw 2 -> tubule, raw 3 -> sheet
                    raw = raw.cpu().numpy()
                    L = np.zeros(raw.shape, dtype=np.uint8)
                    L[raw == 2] = 1   # tubule
                    L[raw == 3] = 2   # sheet
                    L[raw == 1] = 3   # sheet-based tubule
                    _result[0] = (R.cpu().numpy().astype(np.float32), L)

                except Exception:
                    _result[1] = traceback.format_exc()

            import threading
            t0 = time.time()
            th = threading.Thread(target=_run, daemon=True)
            th.start()
            th.join(timeout=timeout)

            if th.is_alive():
                raise RuntimeError(f'inference timed out after {timeout:.0f} s')
            if _result[1] is not None:
                raise RuntimeError(_result[1])

            R, L = _result[0]
            sio.savemat(str(res_file), {'R': R, 'L': L}, format='5')
            print(f'[server] done {time.time()-t0:.2f}s', flush=True)

        except Exception:
            err_file.write_text(traceback.format_exc())
            print(f'[server] ERROR {req_file.name}:\n{traceback.format_exc()}',
                  flush=True)
        finally:
            try:
                req_file.unlink()
            except Exception:
                pass

    time.sleep(0.05)   # 50 ms poll interval

pid_file.unlink(missing_ok=True)
