classdef TestVesselness < matlab.unittest.TestCase
% TestVesselness  Unit tests for vesselnessFilter2DOptimised.
%
% NOTE: FrangiFilter2D and vesselnessFilter2DOptimised use different option
% schemas (BlackWhite vs WhiteOnDark, different Hessian backends) and are
% NOT numerically identical.  Tests here validate the optimised function on
% its own — correctness of output size, range, and line-detection behaviour.
%
% USAGE
%   results = runtests('Vesselness/tests/TestVesselness');
%   table(results)

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            vRoot       = fullfile(fileparts(mfilename('fullpath')), '..');
            projectRoot = fullfile(vRoot, '..');
            addpath(fullfile(vRoot, 'src'));
            addpath(vRoot);                                         % +vesselness package parent
            addpath(fullfile(projectRoot, 'Hessian', 'src'));      % applyHessian2D, eig2image
            addpath(fullfile(projectRoot, 'Helper functions'));     % makeReticulateTestImage
            tc.assumeTrue( ...
                exist('vesselnessFilter2DOptimised', 'file') == 2, ...
                'vesselnessFilter2DOptimised not found on path.');
        end
    end

    % =====================================================================
    % Shared options helper
    % =====================================================================
    methods (Static, Access = private)
        function opts = baseOpts(scaleRange)
            if nargin < 1, scaleRange = [1 3]; end
            opts.FrangiScaleRange = scaleRange;
            opts.FrangiScaleRatio = 1;
            opts.FrangiBetaOne    = 0.5;
            opts.FrangiBetaTwo    = 15;
            opts.WhiteOnDark      = true;   % bright vessels on dark background
            opts.verbose          = false;
        end

        function I = lineImage()
            % Bright horizontal line (sigma 2 px, width ~6 px)
            [~, y] = meshgrid(1:128, 1:128);
            I = single(exp(-(y - 64).^2 / (2*2^2)));
            I = I / max(I(:));
        end
    end

    % =====================================================================
    % Tests
    % =====================================================================
    methods (Test)

        function testOutputSize(tc)
            I = makeReticulateTestImage(128, 128, 42);
            [outOpt, scaleOpt, dirOpt] = vesselnessFilter2DOptimised( ...
                double(I), TestVesselness.baseOpts());
            tc.verifySize(outOpt,   size(I));
            tc.verifySize(scaleOpt, size(I));
            tc.verifySize(dirOpt,   size(I));
        end

        function testResponseRange(tc)
            % Vesselness output must lie in [0, 1]
            I = makeReticulateTestImage(128, 128, 7);
            [outOpt, ~, ~] = vesselnessFilter2DOptimised( ...
                double(I), TestVesselness.baseOpts());
            tc.verifyGreaterThanOrEqual(min(outOpt(:)), 0);
            tc.verifyLessThanOrEqual(   max(outOpt(:)), 1 + 1e-6);
        end

        function testLineDetection(tc)
            % Vesselness must respond more strongly at a bright horizontal
            % line than in the flat background corners
            I    = TestVesselness.lineImage();
            opts = TestVesselness.baseOpts([1 4]);
            [V, ~, ~] = vesselnessFilter2DOptimised(double(I), opts);
            lineResp = mean(V(60:68, 20:108), 'all');
            bgResp   = mean(V(1:20,  1:20),   'all');
            tc.verifyGreaterThan(lineResp, 5 * bgResp + 1e-6);
        end

        function testDeterminism(tc)
            % Identical call must produce identical output
            I = makeReticulateTestImage(128, 128, 99);
            opts = TestVesselness.baseOpts();
            [A, ~, ~] = vesselnessFilter2DOptimised(double(I), opts);
            [B, ~, ~] = vesselnessFilter2DOptimised(double(I), opts);
            tc.verifyEqual(A, B, 'AbsTol', 0);
        end

    end

end
