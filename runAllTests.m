% runAllTests.m
% Run the full test suite for all CurvilinearFilters modules.
%
% All test classes are matlab.unittest.TestCase subclasses and are
% fully discoverable by the MATLAB Test Browser (R2019b+).
%
% Usage (from project root):
%   results = runAllTests();
%
% Returns a struct array of matlab.unittest.TestResult, one per suite.
% A summary table is printed to the command window.

function allResults = runAllTests()

root = fileparts(mfilename('fullpath'));

% ---------------------------------------------------------------------------
% Path setup — add everything needed by all test suites.
% (Each TestCase class also adds its own paths via TestClassSetup, so
%  these are only needed when running runAllTests directly rather than
%  invoking a class independently.)
% ---------------------------------------------------------------------------
addpath(fullfile(root, 'Helper functions'));

% Hessian
addpath(fullfile(root, 'Hessian', 'src'));
addpath(fullfile(root, 'Hessian', 'src', 'engine'));
addpath(fullfile(root, 'Hessian', 'src', 'utils'));
addpath(fullfile(root, 'Hessian', 'src', 'neuriteness'));
addpath(fullfile(root, 'Hessian', 'tests', 'utilities'));


% MFAT
addpath(fullfile(root, 'MFAT'));                       % parent of +mfat package
addpath(fullfile(root, 'MFAT', 'src', 'drivers'));
addpath(fullfile(root, 'MFAT', 'src', 'core'));
addpath(fullfile(root, 'MFAT', 'src', 'responses'));
addpath(fullfile(root, 'MFAT', 'src', 'modifiers'));
addpath(fullfile(root, 'MFAT', 'src', 'utils'));
addpath(fullfile(root, 'MFAT', 'src', 'original'));
addpath(fullfile(root, 'MFAT', 'config'));

% PhaseCongruency
addpath(fullfile(root, 'PhaseCongruency', 'src'));
addpath(fullfile(root, 'PhaseCongruency', 'utils'));

% SOAGK
addpath(fullfile(root, 'SOAGK', 'src'));
addpath(fullfile(root, 'SOAGK', 'utils'));

% Vesselness
addpath(fullfile(root, 'Vesselness', 'src'));
addpath(fullfile(root, 'Vesselness'));         % +vesselness package parent

% BowlerHat
addpath(fullfile(root, 'BowlerHat', 'src'));

% nERdy
addpath(fullfile(root, 'nERdy', 'src'));

% ---------------------------------------------------------------------------
% Test suites — all are now matlab.unittest.TestCase classes
% ---------------------------------------------------------------------------
suites = {
    'Hessian',         fullfile(root, 'Hessian',         'tests', 'TestHessian.m');
    'MFAT',            fullfile(root, 'MFAT',            'tests', 'TestMfat.m');
    'PhaseCongruency', fullfile(root, 'PhaseCongruency', 'tests', 'TestPhaseCongruency.m');
    'SOAGK',           fullfile(root, 'SOAGK',           'tests', 'TestSoagk.m');
    'Vesselness',      fullfile(root, 'Vesselness',      'tests', 'TestVesselness.m');
    'BowlerHat',       fullfile(root, 'BowlerHat',       'tests', 'TestBowlerHat.m');
    'nERdy',           fullfile(root, 'nERdy',           'tests', 'TestNERdy.m');
};

fprintf('\n');
fprintf('%-20s  %6s  %6s  %6s\n', 'Module', 'Pass', 'Fail', 'Total');
fprintf('%s\n', repmat('-', 1, 44));

allResults = [];
anyFailed  = false;

set(groot, 'DefaultFigureVisible', 'off');   % suppress any accidental figure windows

for k = 1:size(suites, 1)
    name = suites{k, 1};
    src  = suites{k, 2};

    try
        r = runtests(src, 'Verbosity', 0);
    catch ME
        fprintf('%-20s  ERROR: %s\n', name, ME.message);
        anyFailed = true;
        continue;
    end

    nPass = sum([r.Passed]);
    nFail = sum([r.Failed]) + sum([r.Incomplete]);
    nTot  = numel(r);

    if nFail > 0
        status    = ' <-- FAILED';
        anyFailed = true;
    else
        status = '';
    end

    fprintf('%-20s  %6d  %6d  %6d%s\n', name, nPass, nFail, nTot, status);
    allResults = [allResults, r]; %#ok<AGROW>
end

set(groot, 'DefaultFigureVisible', 'on');

fprintf('%s\n', repmat('-', 1, 44));
if anyFailed
    fprintf('RESULT: some tests FAILED\n\n');
else
    fprintf('RESULT: all tests passed\n\n');
end

if nargout == 0
    clear allResults;
end
end
