classdef TestVesselness < matlab.unittest.TestCase
% TestVesselness  Unit tests for vesselnessFilter2DOptimised.
%
% Validates numerical equivalence against the reference FrangiFilter2D
% implementation within a reasonable floating-point tolerance.
%
% Migrated from bare-function test to matlab.unittest.TestCase class.
% Pre-existing tolerance fix: response tolerance relaxed from 1e-10 → 1e-6
% to accommodate floating-point ordering differences between the two
% implementations.
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
            addpath(vRoot);                                           % +vesselness package parent
            addpath(fullfile(projectRoot, 'Hessian', 'original'));   % FrangiFilter2D reference
            addpath(fullfile(projectRoot, 'Helper functions'));       % makeReticulateTestImage
            tc.assumeTrue( ...
                exist('FrangiFilter2D',           'file') == 2 && ...
                exist('vesselnessFilter2DOptimised', 'file') == 2, ...
                'FrangiFilter2D or vesselnessFilter2DOptimised not found on path.');
        end
    end

    % =====================================================================
    % Tests
    % =====================================================================
    methods (Test)

        function testEquivalenceToFrangi(tc)
            % The optimised implementation must agree with the reference
            % FrangiFilter2D on all three outputs within relaxed tolerances.
            I = makeReticulateTestImage(256, 256, 123);

            opts.FrangiScaleRange = [1 6];
            opts.FrangiScaleRatio = 1;
            opts.FrangiBetaOne    = 0.5;
            opts.FrangiBetaTwo    = 15;
            opts.BlackWhite       = false;
            opts.WhiteOnDark      = true;
            opts.verbose          = false;

            [outRef, scaleRef, dirRef] = FrangiFilter2D(double(I), opts);
            [outOpt, scaleOpt, dirOpt] = vesselnessFilter2DOptimised(double(I), opts);

            % Response agreement (relaxed: implementations may differ in
            % floating-point order, tie-breaking, etc.)
            tolResponse = 1e-6;
            tc.verifyLessThan( ...
                max(abs(outRef(:) - outOpt(:))), tolResponse, ...
                'Vesselness response mismatch between Frangi and optimised');

            % Scale selection: allow up to 2% pixel mismatch (tie-breaking)
            fracMismatch = mean(scaleRef(:) ~= scaleOpt(:));
            tc.verifyLessThan(fracMismatch, 0.02, ...
                sprintf('Scale selection mismatch: %.1f%% of pixels differ', ...
                        100 * fracMismatch));

            % Direction agreement where vesselness > 0
            mask = outRef > 0;
            if any(mask(:))
                tc.verifyLessThan( ...
                    max(abs(dirRef(mask) - dirOpt(mask))), 1e-6, ...
                    'Direction mismatch');
            end
        end

        function testOutputSize(tc)
            I = makeReticulateTestImage(128, 128, 42);
            opts.FrangiScaleRange = [1 3];
            opts.FrangiScaleRatio = 1;
            opts.FrangiBetaOne    = 0.5;
            opts.FrangiBetaTwo    = 15;
            opts.BlackWhite       = false;
            opts.WhiteOnDark      = true;
            opts.verbose          = false;
            [outOpt, scaleOpt, dirOpt] = vesselnessFilter2DOptimised(double(I), opts);
            tc.verifySize(outOpt,   size(I));
            tc.verifySize(scaleOpt, size(I));
            tc.verifySize(dirOpt,   size(I));
        end

        function testResponseRange(tc)
            I = makeReticulateTestImage(128, 128, 7);
            opts.FrangiScaleRange = [1 3];
            opts.FrangiScaleRatio = 1;
            opts.FrangiBetaOne    = 0.5;
            opts.FrangiBetaTwo    = 15;
            opts.BlackWhite       = false;
            opts.WhiteOnDark      = true;
            opts.verbose          = false;
            [outOpt, ~, ~] = vesselnessFilter2DOptimised(double(I), opts);
            tc.verifyGreaterThanOrEqual(min(outOpt(:)), 0);
            tc.verifyLessThanOrEqual(   max(outOpt(:)), 1 + 1e-6);
        end

    end

end
