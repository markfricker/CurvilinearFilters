classdef TestPhaseCongruency < matlab.unittest.TestCase
% TestPhaseCongruency  Unit tests for phaseCong3Optimised.
%
% Migrated from functiontests pattern to matlab.unittest.TestCase class.
%
% USAGE
%   results = runtests('PhaseCongruency/tests/TestPhaseCongruency');
%   table(results)

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            pcRoot      = fullfile(fileparts(mfilename('fullpath')), '..');
            projectRoot = fullfile(pcRoot, '..');
            addpath(fullfile(pcRoot, 'src'));
            addpath(fullfile(pcRoot, 'utils'));
            addpath(fullfile(projectRoot, 'Helper functions'));
            tc.assumeTrue( ...
                exist('phaseCong3Optimised', 'file') == 2, ...
                'phaseCong3Optimised.m not found on the MATLAB path.');
        end
    end

    % =====================================================================
    % Tests
    % =====================================================================
    methods (Test)

        function testBasicOutputsAndRanges(tc)
            im  = makeReticulateTestImage(256, 256, 123);
            tol = 1e-8;
            [M, m, ori, ft, PC, EO, T, pcSum, noiseMask] = phaseCong3Optimised( ...
                im, 'precision', 'double', 'storeEO', true, 'storePC', true);

            % Sizes
            tc.verifySize(M,         size(im));
            tc.verifySize(m,         size(im));
            tc.verifySize(ori,       size(im));
            tc.verifySize(ft,        size(im));
            tc.verifySize(pcSum,     size(im));
            tc.verifySize(noiseMask, size(im));

            % Finiteness
            tc.verifyFalse(any(~isfinite(M(:))),     'M contains NaN/Inf');
            tc.verifyFalse(any(~isfinite(m(:))),     'm contains NaN/Inf');
            tc.verifyFalse(any(~isfinite(ft(:))),    'featType contains NaN/Inf');
            tc.verifyFalse(any(~isfinite(pcSum(:))), 'pcSum contains NaN/Inf');

            % EO / PC shape
            tc.verifyTrue(iscell(EO));
            tc.verifyEqual(size(EO, 2), 6);   % default norient=6
            tc.verifyEqual(size(EO, 1), 4);   % default nscale=4
            tc.verifyTrue(iscell(PC));
            tc.verifyEqual(numel(PC), 6);

            % Value ranges
            tc.verifyGreaterThanOrEqual(M, -tol);
            tc.verifyGreaterThanOrEqual(m, -tol);
            tc.verifyGreaterThanOrEqual(M - m, -tol);
            tc.verifyGreaterThanOrEqual(ori, 0);
            tc.verifyLessThanOrEqual(   ori, 180);

            for o = 1:numel(PC)
                tc.verifyFalse(any(~isfinite(PC{o}(:))), ...
                    sprintf('PC{%d} contains NaN/Inf', o));
                tc.verifyGreaterThanOrEqual(PC{o}, -1e-6);
                tc.verifyLessThanOrEqual(   PC{o}, 1 + 1e-6);
            end

            tc.verifyTrue(islogical(noiseMask), 'noiseMask should be logical');
            tc.verifyTrue(isscalar(T),          'T should be a scalar');
        end

        function testDefaultStoragePolicy(tc)
            im = makeReticulateTestImage(128, 128, 1);

            % storePC=false by default → PC should be empty
            [~, ~, ~, ~, PC, EO] = phaseCong3Optimised(im, 'precision', 'double');
            tc.verifyTrue(isempty(PC), 'PC should be empty by default (storePC=false)');
            tc.verifyTrue(iscell(EO),  'EO should be a cell when requested as output');

            % storeEO=false → EO empty
            [~, ~, ~, ~, ~, EO2] = phaseCong3Optimised(im, 'precision', 'double', 'storeEO', false);
            tc.verifyTrue(isempty(EO2), 'EO should be [] when storeEO=false');

            % storePC=true → PC is a cell
            [~, ~, ~, ~, PC3] = phaseCong3Optimised(im, 'precision', 'double', 'storePC', true);
            tc.verifyTrue(iscell(PC3), 'PC should be a cell when storePC=true');
        end

        function testDeterminism(tc)
            im = makeReticulateTestImage(192, 192, 77);
            A  = phaseCong3Optimised(im, 'precision', 'double', 'storeEO', false, 'storePC', false);
            B  = phaseCong3Optimised(im, 'precision', 'double', 'storeEO', false, 'storePC', false);
            tc.verifyEqual(A, B, 'AbsTol', 0, ...
                'Outputs should be bitwise identical for identical inputs');
        end

        function testSingleVsDoubleClose(tc)
            imD = makeReticulateTestImage(192, 192, 9);
            imS = single(imD);
            [M1, m1, or1, ft1, ~, ~, T1, pcSum1] = phaseCong3Optimised( ...
                imD, 'precision', 'double', 'storeEO', false, 'storePC', false);
            [M2, m2, or2, ft2, ~, ~, T2, pcSum2] = phaseCong3Optimised( ...
                imS, 'precision', 'single', 'storeEO', false, 'storePC', false);

            tc.verifyEqual(M1,     double(M2),     'RelTol', 5e-3, 'AbsTol', 5e-4);
            tc.verifyEqual(m1,     double(m2),     'RelTol', 5e-3, 'AbsTol', 5e-4);
            tc.verifyEqual(or1,    double(or2),    'AbsTol', 1);
            tc.verifyEqual(ft1,    double(ft2),    'RelTol', 5e-3, 'AbsTol', 5e-4);
            tc.verifyEqual(pcSum1, double(pcSum2), 'RelTol', 5e-3, 'AbsTol', 5e-4);
            tc.verifyEqual(T1,     double(T2),     'RelTol', 5e-3, 'AbsTol', 5e-4);
        end

        function testNoiseMaskSanity(tc)
            im = makeReticulateTestImage(128, 128, 101);
            [~, ~, ~, ~, ~, ~, T, ~, noiseMask] = phaseCong3Optimised( ...
                im, 'precision', 'double', 'noiseMethod', 0.2, ...
                'storeEO', false, 'storePC', false);
            tc.verifyTrue(isscalar(T));
            tc.verifyTrue(islogical(noiseMask));
            tc.verifyGreaterThanOrEqual(nnz(noiseMask), 0);
            tc.verifyLessThanOrEqual(   nnz(noiseMask), numel(noiseMask));

            [~, ~, ~, ~, ~, ~, ~, ~, noiseMaskLow] = phaseCong3Optimised( ...
                im, 'precision', 'double', 'noiseMethod', 0.0, ...
                'storeEO', false, 'storePC', false);
            tc.verifyLessThanOrEqual(nnz(noiseMaskLow), nnz(noiseMask) + 1);
        end

        function testPersistentCacheLoop(tc)
            im  = makeReticulateTestImage(256, 256, 202);
            ref = phaseCong3Optimised(im, 'precision', 'double', 'storeEO', false, 'storePC', false);
            for i = 1:10 %#ok<NASGU>
                out = phaseCong3Optimised(im, 'precision', 'double', ...
                    'storeEO', false, 'storePC', false);
                tc.verifyEqual(out, ref, 'AbsTol', 0);
            end
            clear phaseCong3Optimised
            out2 = phaseCong3Optimised(im, 'precision', 'double', ...
                'storeEO', false, 'storePC', false);
            tc.verifyEqual(out2, ref, 'AbsTol', 0);
        end

    end

end
