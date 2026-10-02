function [R, L] = ernetEnhance(I, params)
% ernetEnhance  Deep-learning ER segmenter via ERnet (vision transformer)
%
% USAGE
%   R = ernetEnhance(I)
%   R = ernetEnhance(I, params)
%   [R, L] = ernetEnhance(...)      % also ERnet's 4-class map
%
% INPUTS
%   I       - 2-D grayscale image (any numeric class); converted to
%             single [0,1] internally.
%   params  - optional struct with fields:
%            .pythonExe   = ''          % Full path to Python executable.
%                                       % '' = auto-detect: checks
%                                       %   %USERPROFILE%\venvs\nerdy first
%                                       %   (ERnet shares the nERdy+ venv),
%                                       %   then pyenv(), then PATH.
%            .device      = 'auto'      % 'auto', 'cuda', or 'cpu'.
%            .threshold   = NaN         % NaN = most likely class (argmax).
%                                       % A value in [0,1] is a fixed
%                                       % threshold on P(ER) = 1 - P(background).
%            .normalize   = true        % API consistency; R is always binary.
%            .useCache    = true        % reuse the result for an identical
%                                       % image + threshold (see CACHE)
%
% OUTPUTS
%   R       - binary enhancement map, single precision, same size as I.
%             Pixels in any ER class (tubule, sheet, sheet-based tubule)
%             = 1; background = 0. Compatible with the [0,1] convention of
%             the other enhancers and usable directly as a segmentation mask.
%   L       - uint8 class map, same size as I: 0 background, 1 tubule,
%             2 sheet, 3 sheet-based tubule (ERnet's most likely class,
%             independent of params.threshold). Scale-sensitive -- see
%             OVERVIEW before using the tubule/sheet split.
%
% CACHE
%   The cisternae step ('ERnet' method) and the Enhance step run ERnet on
%   the same frames, so the last 64 results are kept, keyed by an MD5 of
%   the image bytes, size and threshold (so a changed image can never hit
%   a stale entry). 'clear ernetEnhance' empties it.
%
% REQUIREMENTS
%   The nERdy+ venv at %USERPROFILE%\venvs\nerdy (torch, scipy, numpy) plus
%   einops, and the ERnet model files in %USERPROFILE%\venvs\nerdy\ernet.
%   setupPythonEnvs('ernet') installs einops and downloads the model files
%   from the ERnet-v2 GitHub repository (GPL-3.0, so they are not bundled).
%
% OVERVIEW
%   Uses a persistent ERnet server (../python/ernetServer.py) that loads
%   the model once per MATLAB session, then handles all requests via temp
%   files. The server is started automatically on the first call via the
%   Windows Task Scheduler (same mechanism as nERdyEnhance/cellposeSegment).
%
%   The network is ERnet's single-frame SwinIR model, the only ERnet model
%   with released weights; the temporal (multi-frame) model in the paper is
%   not available. Its classes are background, tubule, sheet and
%   sheet-based tubule; only ER vs background is returned, because the
%   tubule/sheet split changes strongly with image scale (Airyscan ER at
%   0.5x-1.5x: tubule share of ER 69% -> 12%) while the ER/background split
%   does not (11.8-13.3% of the image).
%
% NOTES
%   - For 3-D stacks, call funcEnhanceRun which iterates over Z/T slices.
%   - ERnet does not reject non-ER signal: on a mitochondrial channel it
%     labels the puncta as ER. Use it on an ER channel only.
%   - CPU inference takes ~7 s for 488x488 and ~40-75 s for 1024x1024;
%     params.device = 'auto' uses the GPU when torch sees one.
%   - Check %TEMP%\ernet_work\server_startup.log if the server fails to start.
%
% REFERENCES
%   Lu M, Christensen CN, Weber JM, Konno T, Laubli NF, Scherer KM,
%   Avezov E, Lio P, Lapkin AA, Kaminski Schierle GS, Kaminski CF (2023)
%   ERnet: a tool for the semantic segmentation and quantitative analysis
%   of endoplasmic reticulum topology. Nature Methods 20:569-579.
%   https://doi.org/10.1038/s41592-023-01815-0
%
%   GitHub: https://github.com/charlesnchr/ERnet-v2
%
% EXAMPLE
%   R = ernetEnhance(I);
%
%   p.device    = 'cpu';
%   p.threshold = 0.5;
%   R = ernetEnhance(I, p);
%
% See also: nERdyEnhance, setupPythonEnvs

% --- defaults ---------------------------------------------------------------
if nargin < 2, params = struct(); end
if ~isfield(params, 'pythonExe'),  params.pythonExe  = '';     end
if ~isfield(params, 'device'),     params.device     = 'auto'; end
if ~isfield(params, 'threshold'),  params.threshold  = NaN;   end
if ~isfield(params, 'useCache'),   params.useCache   = true;  end

% --- input validation -------------------------------------------------------
if ~ispc
    error('ernetEnhance:unsupportedPlatform', ...
          ['ernetEnhance: the ERnet server is currently supported on Windows only ' ...
           '(it is launched via Task Scheduler from a Windows venv).']);
end
if size(I, 3) > 1
    error('ernetEnhance:badInput', ...
          'ernetEnhance: expected 2-D grayscale image, got %d-channel input.', ...
          size(I, 3));
end

I = im2single(I);

% --- cache lookup -------------------------------------------------------------
persistent cache cacheKeys
if isempty(cache), cache = containers.Map(); cacheKeys = {}; end
key = '';
if params.useCache
    key = ernetCacheKey(I, params.threshold);
    if isKey(cache, key)
        hit = cache(key);
        R = hit.R;  L = hit.L;
        return;
    end
end

% --- run via persistent server ----------------------------------------------
[R, L] = ernetRunViaServer(I, params);

if params.useCache
    cache(key) = struct('R', R, 'L', L);
    cacheKeys{end+1} = key;
    if numel(cacheKeys) > 64
        remove(cache, cacheKeys{1});
        cacheKeys(1) = [];
    end
end
end


function key = ernetCacheKey(I, threshold)
% keyHash is built in (no Java, which R2026b may lack). It is 64-bit and only
% stable within a session, which is all this in-memory cache needs.
key = sprintf('%016x', keyHash({class(I), size(I), threshold, I}));
end


% ============================================================================
% Server communication
% ============================================================================

function [R, L] = ernetRunViaServer(I, params)

serverScript = ernetResolveServerScript();
workDir      = fullfile(tempdir, 'ernet_work');
if ~exist(workDir, 'dir'), mkdir(workDir); end

pidFile   = fullfile(workDir, 'server.pid');
readyFile = fullfile(workDir, 'server.ready');

if ~ernetServerAlive(pidFile)
    ernetStartServer(params, serverScript, workDir, pidFile, readyFile);
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
        error('ernetEnhance:timeout', ...
              'ERnet server timed out after %d s.', pollTimeout);
    end
    pause(0.25);
end

if exist(errFile, 'file')
    msg = fileread(errFile);
    delete(errFile);
    error('ernetEnhance:serverError', 'ERnet server error:\n%s', msg);
end

result = load(resFile);
delete(resFile);
R = single(result.R);
if isfield(result, 'L')
    L = uint8(result.L);
    return;
end
% A server started by an older ernetServer.py (it outlives MATLAB) returns
% no class map: ask it to exit, wait for the process to end, and re-run --
% the next call starts the current script.
if isfield(params, 'restarted') && params.restarted
    error('ernetEnhance:staleServer', ...
          'ERnet server still returns no class map after a restart. Check %s.', ...
          fullfile(workDir, 'server_startup.log'));
end
fid = fopen(fullfile(workDir, 'exit.req'), 'w'); fclose(fid);
t0 = tic;
while ernetServerAlive(pidFile) && toc(t0) < 15
    pause(0.25);
end
if ernetServerAlive(pidFile)
    pid = str2double(strtrim(fileread(pidFile)));
    system(sprintf('taskkill /PID %d /F > NUL 2>&1', pid));
end
params.restarted = true;
[R, L] = ernetRunViaServer(I, params);
end


function ernetStartServer(params, serverScript, workDir, pidFile, readyFile)
logFile   = fullfile(workDir, 'server_startup.log');
errorFile = fullfile(workDir, 'server.error');

% Remove stale files from any previous run
for f = {readyFile, pidFile, errorFile}
    if exist(f{1}, 'file'), delete(f{1}); end
end

pyExe = ernetFindPython(params);

taskName = 'MATLABERnetServer';
tr = sprintf('"%s" "%s" "%s"', pyExe, serverScript, workDir);
system(sprintf('schtasks /Delete /TN "%s" /F > NUL 2>&1', taskName));
createCmd = sprintf( ...
    'schtasks /Create /F /TN "%s" /TR "%s" /SC ONCE /SD 01/01/2000 /ST 00:00', ...
    taskName, strrep(tr, '"', '\"'));
system(createCmd);
% schtasks /Create defaults DisallowStartIfOnBatteries = TRUE which silently
% parks the task as "Queued" on a laptop running on battery.  Flip both flags
% off before running, matching the fix in cellposeSegment.m.
powerCmd = sprintf(['powershell -NoProfile -Command ' ...
    '"Set-ScheduledTask -TaskName ''%s'' -Settings ' ...
    '(New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries ' ...
    '-DontStopIfGoingOnBatteries -MultipleInstances IgnoreNew)"'], taskName);
system(powerCmd);
system(sprintf('schtasks /Run /TN "%s"', taskName));

% --- waitbar with Cancel ------------------------------------------------------
wb = waitbar(0, 'Starting ERnet server...', 'Name', 'ERnet', ...
             'CreateCancelBtn', @(~,~) setappdata(gcbf, 'cancel', true));
setappdata(wb, 'cancel', false);
wbClean = onCleanup(@() ernetCloseWaitbar(wb));

% Wait for PID file (up to 60 s — covers slow Task Scheduler launch)
t0 = tic;
while ~ernetServerAlive(pidFile) && toc(t0) < 60
    if ~ishandle(wb) || getappdata(wb, 'cancel')
        error('ernetEnhance:cancelled', 'ERnet server startup cancelled.');
    end
    if exist(errorFile, 'file')
        msg = fileread(errorFile);
        error('ernetEnhance:serverCrash', ...
              'ERnet server crashed during startup:\n%s', msg);
    end
    waitbar(min(toc(t0)/60, 0.4), wb, 'Starting ERnet server...');
    pause(0.5);
end
if ~ernetServerAlive(pidFile)
    error('ernetEnhance:serverTimeout', ...
          ['ERnet server did not start within 60 s.\n' ...
           'Check log: %s\n' ...
           'Or start manually from PowerShell:\n' ...
           '  python "%s" "%s"'], logFile, serverScript, workDir);
end

% Wait for ready file (model loaded — up to 120 s for CUDA torch)
t0 = tic;
while ~exist(readyFile, 'file') && toc(t0) < 120
    if ~ishandle(wb) || getappdata(wb, 'cancel')
        error('ernetEnhance:cancelled', 'ERnet server startup cancelled.');
    end
    if exist(errorFile, 'file')
        msg = fileread(errorFile);
        error('ernetEnhance:serverCrash', ...
              'ERnet server crashed during startup:\n%s', msg);
    end
    if ~ernetServerAlive(pidFile)
        msg = '';
        if exist(logFile, 'file'), msg = fileread(logFile); end
        error('ernetEnhance:serverCrash', ...
              ['ERnet server process died before model loaded.\n\n' ...
               'Startup log:\n%s'], msg);
    end
    waitbar(0.4 + min(toc(t0)/120, 0.55), wb, 'Loading ERnet model...');
    pause(0.5);
end
if ~exist(readyFile, 'file')
    error('ernetEnhance:modelTimeout', ...
          ['ERnet model did not load within 120 s.\n' ...
           'Check log: %s'], logFile);
end

waitbar(1, wb, 'ERnet server ready.');
pause(0.3);
if ishandle(wb), delete(wb); end
end


function ernetCloseWaitbar(wb)
if ishandle(wb), delete(wb); end
end


function serverScript = ernetResolveServerScript()
% ernetServer.py sits in ../python relative to this file, in the dev layout
% and in the toolbox (buildToolbox stages ERnet/src and ERnet/python).
scriptDir    = fileparts(mfilename('fullpath'));   % .../ERnet/src
serverScript = fullfile(fileparts(scriptDir), 'python', 'ernetServer.py');
if ~isfile(serverScript)
    error('ernetEnhance:notFound', ...
          ['ernetServer.py not found.\n' ...
           'Expected at: %s\n' ...
           'Run setupPythonEnvs(''ernet'') to configure the environment.'], ...
          serverScript);
end
end


function p = ernetFindPython(params)
% Search order: explicit param → nerdy venv (shared) → pyenv → PATH
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


function alive = ernetServerAlive(pidFile)
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
