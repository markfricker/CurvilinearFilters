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
%                                       % '' = use 'python' on system PATH.
%            .nerdyInfer  = ''          % Full path to nerdy_infer.py.
%                                       % '' = auto-detect relative to this
%                                       %      function file.
%            .modelPath   = ''          % Full path to .pth weights file.
%                                       % '' = auto-detect (default weights).
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
%   Python 3.9+ with the following packages installed:
%     torch>=2.0.0, torchvision, scikit-image, pillow, numpy
%   The nERdy+ source tree (with pre-trained weights) must be available on
%   disk.  By default the function looks for:
%     <this_file>/../../../third party/nERdy integration/nerdy_infer.py
%   and for the model weights at:
%     <nerdy_infer_dir>/nERdy/nERdy+/NNet_groupy_p4m_v2_VecAdam.pth
%
% OVERVIEW
%   Segmentation uses the D4-equivariant encoder-decoder (D4nERdy) from the
%   nERdy+ paper (Samudre et al., 2024, Nature Communications Biology).
%   The network was trained on confocal/STED ER fluorescence images and
%   outperforms classical filters for tubule detection.
%
%   The function writes the input slice to a temporary TIFF file, calls
%   nerdy_infer.py as a subprocess (avoiding MATLAB pyenv dependency), reads
%   the resulting binary mask, and returns it as a single-precision [0,1]
%   image.
%
% NOTES
%   - For 3-D stacks, call funcEnhanceRun which iterates over Z/T slices.
%   - GPU inference (params.device = 'cuda') is ~5× faster than CPU for
%     512×512 images on a modern laptop GPU.
%   - The first call per session incurs model-loading overhead (~2 s GPU,
%     ~5 s CPU for the 1.5 M-parameter D4nERdy network).
%   - If the subprocess fails, the function throws an informative error
%     that includes the Python stderr output for easy debugging.
%
% REFERENCES
%   Samudre A, Gao G, Cardoen B, Joshi B, Nabi IR, Hamarneh G (2025)
%   nERdy: network analysis of endoplasmic reticulum dynamics.
%   Communications Biology. https://doi.org/10.1038/s42003-025-08892-1
%
%   GitHub: https://github.com/NanoscopyAI/nERdy
%
% EXAMPLE
%   % Basic ER segmentation with GPU auto-selection
%   R = nERdyEnhance(I);
%
%   % Force CPU, custom weights
%   p.device    = 'cpu';
%   p.modelPath = 'C:\mymodels\nERdy_finetuned.pth';
%   R = nERdyEnhance(I, p);
%
%   % Fixed probability threshold instead of Otsu
%   p.threshold = 0.5;
%   R = nERdyEnhance(I, p);
%
% See also: hessian2DFilters, neuriteness2D, cellposeSegment

% --- defaults ---------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'pythonExe'),  params.pythonExe  = '';     end
if ~isfield(params, 'nerdyInfer'), params.nerdyInfer = '';     end
if ~isfield(params, 'modelPath'),  params.modelPath  = '';     end
if ~isfield(params, 'device'),     params.device     = 'auto'; end
if ~isfield(params, 'threshold'),  params.threshold  = NaN;   end

% --- input validation -------------------------------------------------------
if size(I, 3) > 1
    error('nERdyEnhance:badInput', ...
          'nERdyEnhance: expected 2-D grayscale image, got %d-channel input.', ...
          size(I, 3));
end

I = im2single(I);

% --- resolve Python executable ----------------------------------------------
if isempty(params.pythonExe)
    pyExe = nerdyFindPython();
else
    pyExe = params.pythonExe;
end

% --- resolve nerdy_infer.py path -------------------------------------------
if isempty(params.nerdyInfer)
    scriptDir    = fileparts(mfilename('fullpath'));   % CurvilinearFilters_sandbox/nERdy/src
    nerdyInfer   = fullfile(scriptDir, '..', '..', '..', ...
                            'third party', 'nERdy integration', 'nerdy_infer.py');
    nerdyInfer   = char(java.io.File(nerdyInfer).getCanonicalPath());
else
    nerdyInfer = params.nerdyInfer;
end
if ~isfile(nerdyInfer)
    error('nERdyEnhance:notFound', ...
          ['nerdy_infer.py not found at:\n  %s\n' ...
           'Set params.nerdyInfer to the correct path, or ensure the\n' ...
           'nERdy integration folder is adjacent to BlobFilters_sandbox.'], ...
          nerdyInfer);
end

% --- write input to temp TIFF -----------------------------------------------
tmpIn  = [tempname, '_nerdy_in.tif'];
tmpOut = [tempname, '_nerdy_out.tif'];
cleanupObj = onCleanup(@() nerdyCleanTmp(tmpIn, tmpOut));

imwrite(uint8(I * 255), tmpIn);

% --- build command ----------------------------------------------------------
cmd = sprintf('"%s" "%s" --input "%s" --output "%s" --device %s', ...
              pyExe, nerdyInfer, tmpIn, tmpOut, params.device);

if ~isempty(params.modelPath)
    cmd = [cmd, sprintf(' --model "%s"', params.modelPath)];
end

% Threshold: NaN or negative → use default Otsu postprocessing (no --threshold flag).
% A non-negative finite value is an explicit override in [0,1].
if ~isnan(params.threshold) && isfinite(params.threshold) && params.threshold >= 0
    cmd = [cmd, sprintf(' --threshold %.6f', params.threshold)];
end

% --- run Python subprocess --------------------------------------------------
[status, cmdout] = system(cmd);
if status ~= 0
    error('nERdyEnhance:pythonError', ...
          ['nERdy+ inference failed (exit code %d).\n' ...
           'Command: %s\n' ...
           'Output:\n%s\n\n' ...
           'Tip: run the command above in a terminal to diagnose. Ensure\n' ...
           '  torch, torchvision, scikit-image, pillow are installed in\n' ...
           '  the Python environment at: %s'], ...
          status, cmd, cmdout, pyExe);
end

% --- read result ------------------------------------------------------------
if ~isfile(tmpOut)
    error('nERdyEnhance:noOutput', ...
          ['nERdy+ did not produce an output file.\nExpected: %s\n' ...
           'Command output:\n%s'], tmpOut, cmdout);
end

mask = imread(tmpOut);
if ~islogical(mask) && max(mask(:)) > 1
    R = single(mask > 0);   % binarise from uint8 0/255
else
    R = single(mask);
end
end


% ============================================================================
% Local helpers
% ============================================================================

function p = nerdyFindPython()
% nerdyFindPython  Locate a Python executable suitable for nERdy+ inference.
%
% Search order (first found wins):
%   1. Dedicated nERdy venv  C:\Users\<user>\venvs\nerdy\Scripts\python.exe
%      (CUDA torch, scikit-image — created specifically for nERdy+)
%   2. MATLAB's pyenv executable (if configured and non-empty)
%   3. 'python' on system PATH

% 1. Dedicated nerdy venv alongside the cellpose venv
nerdyVenv = fullfile(getenv('USERPROFILE'), 'venvs', 'nerdy', ...
                     'Scripts', 'python.exe');
if isfile(nerdyVenv)
    p = nerdyVenv;
    return;
end

% 2. MATLAB pyenv
try
    pe = pyenv();
    if ~isempty(pe.Executable)
        p = char(pe.Executable);
        return;
    end
catch
end

% 3. System PATH fallback
p = 'python';
end


function nerdyCleanTmp(varargin)
% Delete temporary files, ignoring errors.
for k = 1:nargin
    f = varargin{k};
    if isfile(f)
        try, delete(f); catch, end
    end
end
end
