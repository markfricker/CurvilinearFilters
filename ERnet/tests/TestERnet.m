classdef TestERnet < matlab.unittest.TestCase
% TestERnet  Unit tests for ernetEnhance.
%
% Tests that do NOT require the ERnet Python environment:
%   testERnet_badInput_3D_errors      – 3-D input throws ernetEnhance:badInput
%   testERnet_serverScriptPresent     – ../python/ernetServer.py ships with the wrapper
%
% Tests that DO require ERnet (skipped when absent — run setupPythonEnvs('ernet')):
%   testERnet_outputSizeClassBinary   – size, single class, only 0s and 1s
%   testERnet_oddSize                 – sizes that are not multiples of the
%                                       SwinIR window (4) are padded and cropped back
%   testERnet_segmentsSyntheticNetwork – a synthetic tubule lattice is found as
%                                       ER and flat background is not
%   testERnet_thresholdMonotone       – a lower P(ER) threshold keeps at least
%                                       as many pixels as a higher one
%
% USAGE
%   results = runtests('ERnet/tests/TestERnet');
%   table(results)

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            ernetRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(ernetRoot, 'src'));
        end
    end

    % =====================================================================
    % Helpers
    % =====================================================================
    methods (Static, Access = private)

        function tf = isERnetAvailable()
            % The nERdy venv with einops + both downloaded model files.
            venv = fullfile(getenv('USERPROFILE'), 'venvs', 'nerdy');
            tf = isfile(fullfile(venv, 'Scripts', 'python.exe')) && ...
                 isfile(fullfile(venv, 'ernet', 'swinir_rcab_arch.py')) && ...
                 isfile(fullfile(venv, 'ernet', '20220306_ER_4class_swinir_nch1.pth'));
            if tf
                [rc, ~] = system(sprintf('"%s" -c "import einops"', ...
                    fullfile(venv, 'Scripts', 'python.exe')));
                tf = rc == 0;
            end
        end

        function I = latticeImage(n)
            % Synthetic ER-like image: a noisy polygonal lattice of thin
            % Gaussian-profile tubules (~2 px sigma) on a dim background.
            if nargin < 1, n = 128; end
            rng(1);
            [x, y] = meshgrid(1:n, 1:n);
            I = zeros(n, 'single');
            s = 1.2;
            for c = 12:24:n
                I = max(I, exp(-(x - c - 3*sin(y/15)).^2 / (2*s^2)));
                I = max(I, exp(-(y - c - 3*cos(x/17)).^2 / (2*s^2)));
            end
            I = single(0.05 + 0.9*I + 0.03*randn(n));
            I = min(max(I, 0), 1);
        end

    end

    % =====================================================================
    % Error-condition tests (no Python required)
    % =====================================================================
    methods (Test)

        function testERnet_badInput_3D_errors(tc)
            I3D = repmat(TestERnet.latticeImage(32), [1 1 3]);
            tc.verifyError(@() ernetEnhance(I3D), 'ernetEnhance:badInput');
        end

        function testERnet_serverScriptPresent(tc)
            srcDir = fileparts(which('ernetEnhance'));
            tc.verifyTrue(isfile(fullfile(srcDir, '..', 'python', 'ernetServer.py')), ...
                'ernetServer.py must sit in ../python next to ernetEnhance.m');
        end

    end

    % =====================================================================
    % Functional tests (skipped when ERnet is not installed)
    % =====================================================================
    methods (Test)

        function testERnet_outputSizeClassBinary(tc)
            tc.assumeTrue(TestERnet.isERnetAvailable(), ...
                'ERnet not installed (setupPythonEnvs(''ernet'')) — test skipped');
            I = TestERnet.latticeImage(128);
            R = ernetEnhance(I, struct('device', 'cpu'));
            tc.verifySize(R, size(I));
            tc.verifyClass(R, 'single');
            tc.verifyTrue(all(R(:) == 0 | R(:) == 1), ...
                'Output must be binary (0s and 1s only)');
        end

        function testERnet_oddSize(tc)
            tc.assumeTrue(TestERnet.isERnetAvailable(), ...
                'ERnet not installed — test skipped');
            I = TestERnet.latticeImage(130);
            I = I(1:127, 1:129);
            R = ernetEnhance(I, struct('device', 'cpu'));
            tc.verifySize(R, [127 129]);
        end

        function testERnet_segmentsSyntheticNetwork(tc)
            tc.assumeTrue(TestERnet.isERnetAvailable(), ...
                'ERnet not installed — test skipped');
            I = TestERnet.latticeImage(128);
            R = ernetEnhance(I, struct('device', 'cpu'));
            onTube  = I > 0.6;
            offTube = I < 0.15;
            tc.verifyGreaterThan(mean(R(onTube)), 0.8, ...
                'Most tubule-centre pixels should be labelled ER');
            tc.verifyLessThan(mean(R(offTube)), 0.2, ...
                'Most background pixels should not be labelled ER');
        end

        function testERnet_thresholdMonotone(tc)
            tc.assumeTrue(TestERnet.isERnetAvailable(), ...
                'ERnet not installed — test skipped');
            I = TestERnet.latticeImage(128);
            Rlo = ernetEnhance(I, struct('device', 'cpu', 'threshold', 0.2));
            Rhi = ernetEnhance(I, struct('device', 'cpu', 'threshold', 0.8));
            tc.verifyGreaterThanOrEqual(nnz(Rlo), nnz(Rhi));
            tc.verifyTrue(all(Rhi(:) <= Rlo(:)), ...
                'Every pixel kept at threshold 0.8 must also be kept at 0.2');
        end

    end

end
