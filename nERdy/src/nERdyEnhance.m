function R = nERdyEnhance(I, params)
% nERdyEnhance  Deep-learning ER tubule segmenter via nERdy+
%
% USAGE
%   R = nERdyEnhance(I)
%   R = nERdyEnhance(I, params)
%
% INPUTS
%   I       - 2-D grayscale image (any numeric class); converted to
%             single [0,1] internally.
%   params  - optional struct with fields:
%            .pythonExe   = ''          % Full path to Python executable.
%                                       % '' = auto-detect: checks
%                                       %   %USERPROFILE%\venvs\nerdy first,
%                                       %   then pyenv(), then PATH.
%            .device      = 'auto'      % 'auto', 'cuda', or 'cpu'.
%            .threshold   = NaN         % Fixed threshold in [0,1] for
%                                       % binarisation.  NaN = use default
%                                       % Otsu-based post-processing.
%            .normalize   = true        % API consistency; R is always binary.
%
% OUTPUTS
%   R       - binary enhancement map, single precision, same size as I.
%             Pixels classified as ER tubules = 1; background = 0.
%             Compatible with the [0,1] convention of all other BlobFilters
%             enhancers and can be used directly as a segmentation mask.
%
% REQUIREMENTS
%   A Python venv at %USERPROFILE%\venvs\nerdy with:
%     torch>=2.0.0, torchvision, scikit-image, pillow, numpy
%
%   To create the venv (run once in PowerShell):
%     python -m venv $env:USERPROFILE\venvs\nerdy
%     $env:USERPROFILE\venvs\nerdy\Scripts\pip install torch torchvision scikit-image pillow numpy
%   For CUDA (adjust cu124 to match your driver):
%     $env:USERPROFILE\venvs\nerdy\Scripts\pip install torch torchvision `
%         --index-url https://download.pytorch.org/whl/cu124
%
%   The nERdy+ source tree (model.py + NNet_groupy_p4m_v2_VecAdam.pth) must
%   be in the nERdy/nERdy+/ subfolder alongside nerdyServer.py.
%
% OVERVIEW
%   Uses a persistent nERdy+ server (nerdyServer.py) that loads the model
%   once per MATLAB session, then handles all requests via temp files.
%   The server is started automatically on the first call via the Windows
%   Task Scheduler (same mechanism as cellposeServer.py).
%
%   This avoids the ~2–5 s model-load overhead that the old subprocess
%   approach incurred on every frame.
%
% NOTES
%   - For 3-D stacks, call funcEnhanceRun which iterates over Z/T slices.
%   - GPU inference (params.device = 'cuda') is ~5× faster than CPU for
%     512×512 images.  The server defaults to CPU (safe in Task Scheduler);
%     set params.device = 'cuda' to switch per request.
%   - The server writes a ready file after model load, so MATLAB waits for
%     the full ~5 s startup on the first call, not just the PID file.
%   - Check %TEMP%\nerdy_work\server_startup.log if the server fails to start.
%
% REFERENCES
%   Samudre A, Gao G, Cardoen B, Joshi B, Nabi IR, Hamarneh G (2025)
%   nERdy: network analysis of endoplasmic reticulum dynamics.
%   Communications Biology. https://doi.org/10.1038/s42003-025-08892-1
%
%   GitHub: https://github.com/NanoscopyAI/nERdy
%
% EXAMPLE
%   R = nERdyEnhance(I);
%
%   p.device    = 'cuda';
%   p.threshold = 0.5;
%   R = nERdyEnhance(I, p);
%
% See also: hessian2DFilters, neuriteness2D, cellposeSegment

% --- defaults ---------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'pythonExe'),  params.pythonExe  = '';     end
if ~isfield(params, 'device'),     params.device     = 'auto'; end
if ~isfield(params, 'threshold'),  params.threshold  = NaN;   end

% --- input validation -------------------------------------------------------
if size(I, 3) > 1
    error('nERdyEnhance:badInput', ...
          'nERdyEnhance: expected 2-D grayscale image, got %d-channel input.', ...
          size(I, 3));
end

I = im2single(I);

% --- run via persistent server ----------------------------------------------
R = nerdyRunViaServer(I, params);
end


% ============================================================================
% Server communication
% ============================================================================

function R = nerdyRunViaServer(I, params)

serverScript = nerdyResolveServerScript();
workDir      = fullfile(tempdir, 'nerdy_work');
if ~exist(workDir, 'dir'), mkdir(workDir); end

pidFile   = fullfile(workDir, 'server.pid');
readyFile = fullfile(workDir, 'server.ready');

if ~nerdyServerAlive(pidFile)
    nerdyStartServer(params, serverScript, workDir, pidFile, readyFile);
end

% Write request (atomic: write to .tmp then rename)
reqId   = sprintf('%s_%d', datestr(now,'yyyymmddHHMMSSFFF'), randi(99999));
reqFile = fullfile(workDir, [reqId '.req.mat']);
resFile = fullfile(workDir, [reqId '.res.mat']);
errFile = fullfile(workDir, [reqId '.err']);

device    = params.device;    %#ok<NASGU>
threshold = params.threshold; %#ok<NASGU>
tmpFile   = [reqFile '.tmp'];
save(tmpFile, 'I', 'device', 'threshold', '-v6');
movefile(tmpFile, reqFile);

% Poll for result
pollTimeout = 300;
t0 = tic;
while ~exist(resFile, 'file') && ~exist(errFile, 'file')
    if toc(t0) > pollTimeout
        try, delete(reqFile); catch, end
        error('nERdyEnhance:timeout', ...
              'nERdy+ server timed out after %d s.', pollTimeout);
    end
    pause(0.25);
end

if exist(errFile, 'file')
    msg = fileread(errFile);
    delete(errFile);
    error('nERdyEnhance:serverError', 'nERdy+ server error:\n%s', msg);
end

result = load(resFile);
delete(resFile);
R = single(result.R);
end


function nerdyStartServer(params, serverScript, workDir, pidFile, readyFile)
logFile   = fullfile(workDir, 'server_startup.log');
errorFile = fullfile(workDir, 'server.error');

% Remove stale files from any previous run
for f = {readyFile, pidFile, errorFile}
    if exist(f{1}, 'file'), delete(f{1}); end
end

pyExe = nerdyFindPython(params);

taskName = 'MATLABNerdyServer';
tr = sprintf('"%s" "%s" "%s"', pyExe, serverScript, workDir);
system(sprintf('schtasks /Delete /TN "%s" /F > NUL 2>&1', taskName));
createCmd = sprintf( ...
    'schtasks /Create /F /TN "%s" /TR "%s" /SC ONCE /SD 01/01/2000 /ST 00:00', ...
    taskName, strrep(tr, '"', '\"'));
system(createCmd);
system(sprintf('schtasks /Run /TN "%s"', taskName));

% --- waitbar with Cancel ------------------------------------------------------
wb = waitbar(0, 'Starting nERdy+ server...', 'Name', 'nERdy+', ...
             'CreateCancelBtn', @(~,~) setappdata(gcbf, 'cancel', true));
setappdata(wb, 'cancel', false);
wbClean = onCleanup(@() nerdyCloseWaitbar(wb));

% Wait for PID file (up to 60 s — covers slow Task Scheduler launch)
t0 = tic;
while ~nerdyServerAlive(pidFile) && toc(t0) < 60
    if ~ishandle(wb) || getappdata(wb, 'cancel')
        error('nERdyEnhance:cancelled', 'nERdy+ server startup cancelled.');
    end
    if exist(errorFile, 'file')
        msg = fileread(errorFile);
        error('nERdyEnhance:serverCrash', ...
              'nERdy+ server crashed during startup:\n%s', msg);
    end
    waitbar(min(toc(t0)/60, 0.4), wb, 'Starting nERdy+ server...');
    pause(0.5);
end
if ~nerdyServerAlive(pidFile)
    error('nERdyEnhance:serverTimeout', ...
          ['nERdy+ server did not start within 60 s.\n' ...
           'Check log: %s\n' ...
           'Or start manually from PowerShell:\n' ...
           '  python "%s" "%s"'], logFile, serverScript, workDir);
end

% Wait for ready file (model loaded — up to 120 s for CUDA torch)
t0 = tic;
while ~exist(readyFile, 'file') && toc(t0) < 120
    if ~ishandle(wb) || getappdata(wb, 'cancel')
        error('nERdyEnhance:cancelled', 'nERdy+ server startup cancelled.');
    end
    if exist(errorFile, 'file')
        msg = fileread(errorFile);
        error('nERdyEnhance:serverCrash', ...
              'nERdy+ server crashed during startup:\n%s', msg);
    end
    if ~nerdyServerAlive(pidFile)
        msg = '';
        if exist(logFile, 'file'), msg = fileread(logFile); end
        error('nERdyEnhance:serverCrash', ...
              ['nERdy+ server process died before model loaded.\n\n' ...
               'Startup log:\n%s'], msg);
    end
    waitbar(0.4 + min(toc(t0)/120, 0.55), wb, 'Loading nERdy+ model...');
    pause(0.5);
end
if ~exist(readyFile, 'file')
    error('nERdyEnhance:modelTimeout', ...
          ['nERdy+ model did not load within 120 s.\n' ...
           'Check log: %s'], logFile);
end

waitbar(1, wb, 'nERdy+ server ready.');
pause(0.3);
if ishandle(wb), delete(wb); end
end


function nerdyCloseWaitbar(wb)
if ishandle(wb), delete(wb); end
end



function serverScript = nerdyResolveServerScript()
if isdeployed()
    serverScript = fullfile(ctfroot, 'nERdy_integration', 'nerdyServer.py');
else
    scriptDir    = fileparts(mfilename('fullpath')); % .../nERdy/src
    serverScript = fullfile(scriptDir, '..', '..', '..', ...
                            'third party', 'nERdy integration', 'nerdyServer.py');
    serverScript = char(java.io.File(serverScript).getCanonicalPath());
end
if ~isfile(serverScript)
    error('nERdyEnhance:notFound', ...
          ['nerdyServer.py not found at:\n  %s\n' ...
           'Ensure the nERdy integration folder is adjacent to ' ...
           'CurvilinearFilters_sandbox.'], serverScript);
end
end


function p = nerdyFindPython(params)
% Search order: explicit param → nerdy venv → pyenv → PATH
if isfield(params, 'pythonExe') && ~isempty(params.pythonExe)
    p = params.pythonExe;
    return;
end
nerdyVenv = fullfile(getenv('USERPROFILE'), 'venvs', 'nerdy', ...
                     'Scripts', 'python.exe');
if isfile(nerdyVenv)
    p = nerdyVenv;
    return;
end
try
    pe = pyenv();
    if ~isempty(pe.Executable)
        p = char(pe.Executable);
        return;
    end
catch
end
p = 'python';
end


function alive = nerdyServerAlive(pidFile)
alive = false;
if ~exist(pidFile, 'file'), return; end
try
    pid = str2double(strtrim(fileread(pidFile)));
    if isnan(pid) || pid <= 0, return; end
    [~, out] = system(sprintf('tasklist /FI "PID eq %d" /NH 2>NUL', pid));
    alive = contains(out, num2str(pid));
catch
end
end
