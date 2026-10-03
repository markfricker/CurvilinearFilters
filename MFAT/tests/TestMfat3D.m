classdef TestMfat3D < matlab.unittest.TestCase
% TestMfat3D  Unit tests for mfatLambda3D / mfatProb3D / mfatCore3D.
%
% USAGE
%   results = runtests('MFAT/tests/TestMfat3D');

    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            mfatRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(mfatRoot, 'src', 'drivers'));
            addpath(fullfile(mfatRoot, 'src', 'core'));
            addpath(fullfile(mfatRoot, 'src', 'responses'));
            addpath(fullfile(mfatRoot, 'src', 'utils'));
            addpath(fullfile(mfatRoot, '..', 'Hessian', 'src'));   % applyHessian3DAniso, eig3volume
        end
    end

    methods (Static, Access = private)
        function V = tube(sz, radius, spacing)
            % bright cylinder along dim 1 through the volume centre,
            % radius in physical units
            if nargin < 3, spacing = [1 1 1]; end
            [~, Y, Z] = ndgrid((1:sz(1))*spacing(1), (1:sz(2))*spacing(2), (1:sz(3))*spacing(3));
            c2 = sz(2)/2*spacing(2);  c3 = sz(3)/2*spacing(3);
            V = double((Y-c2).^2 + (Z-c3).^2 <= radius^2);
            V = imgaussfilt3(V, 0.7);
        end
    end

    methods (Test)

        function testOutputRangeAndSize(tc)
            V = TestMfat3D.tube([32 32 32], 2.5);
            out = mfatLambda3D(V, 1:0.5:2.5);
            tc.verifySize(out, size(V));
            tc.verifyClass(out, 'single');
            tc.verifyGreaterThanOrEqual(min(out(:)), 0);
            tc.verifyLessThanOrEqual(max(out(:)), 1 + 1e-6);
            p = mfatProb3D(V, 1:0.5:2.5);
            tc.verifyGreaterThanOrEqual(min(p(:)), 0);
            tc.verifyLessThanOrEqual(max(p(:)), 1 + 1e-6);
        end

        function testTubeVsSheetDiscrimination(tc)
            % the decisive 3D property: a tube outscores a sheet (no 2D
            % per-slice filter can separate these). Both in ONE volume --
            % the output is rescaled by its own max, so a sheet-only volume
            % would be stretched to ~1 and say nothing.
            % MFAT-λ: a sheet's eigenvalues (0,a,a) give sqrt(3/2)*FA~0.71,
            % i.e. a ~0.3 floor (measured: tube 0.84, sheet 0.35); the
            % Beta-LLR in mfatProb3D pushes that floor down (0.76 vs 0.04).
            sz = [48 48 48];
            [~, Y, Z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
            tube  = double((Y-14).^2 + (Z-24).^2 <= 2.5^2);
            sheet = zeros(sz);  sheet(:, 30:48, 22:26) = 1;
            V = imgaussfilt3(tube + sheet, 0.7);
            tCore = imgaussfilt3(tube, 0.7) > 0.5;
            sCore = imgaussfilt3(sheet, 0.7) > 0.5;
            sCore(:, 1:34, :) = false;                    % away from the sheet's edge
            bg = ~imdilate(tCore | sCore, ones(5,5,5));
            sig = 1:0.5:2.5;

            o = mfatLambda3D(V, sig);
            tc.verifyGreaterThan(mean(o(tCore)), 0.6, 'tube core should respond strongly');
            tc.verifyGreaterThan(mean(o(tCore)), 2*mean(o(sCore)), 'MFAT-λ: tube > 2x sheet');
            tc.verifyLessThan(mean(o(bg)), 0.05, 'background should stay low');

            p = mfatProb3D(V, sig);
            tc.verifyGreaterThan(mean(p(tCore)), 5*mean(p(sCore)), 'MFAT-Prob: tube > 5x sheet');
        end

        function testAxialTube_nativeAnisotropicGrid(tc)
            % tube running along the coarse Z axis, processed on the native
            % grid with per-axis physical sigma -- invisible to a per-slice
            % 2D ridge filter, must be found here
            spacing = [0.1 0.1 0.3];
            sz = [40 40 24];
            [X, Y, ~] = ndgrid((1:sz(1))*spacing(1), (1:sz(2))*spacing(2), (1:sz(3))*spacing(3));
            V = double((X-2).^2 + (Y-2).^2 <= 0.25^2);
            V = imgaussfilt3(V, [0.7 0.7 0.3]);
            sig = {[0.15 0.15 0.3], [0.25 0.25 0.3]};
            out = mfatLambda3D(V, sig, 'spacing', spacing);
            core = V > 0.5*max(V(:));
            tc.verifyGreaterThan(mean(out(core)), 0.5);
            tc.verifyGreaterThan(mean(out(core)), 10*mean(out(~imdilate(core, ones(7,7,3)))));
        end

        function testPolarityGate(tc)
            V = TestMfat3D.tube([32 32 32], 2.5);
            Vd = 1 - V;                                   % dark tube on bright background
            sig = 1:0.5:2.5;
            core = V > 0.5;
            oDefault = mfatLambda3D(Vd, sig);
            oFlip    = mfatLambda3D(Vd, sig, 'whiteOnDark', false);
            tc.verifyLessThan(mean(oDefault(core)), 0.1, 'dark tube must be suppressed by default');
            tc.verifyGreaterThan(mean(oFlip(core)), 0.5, 'whiteOnDark=false must recover it');
        end

        function testCandidateMaskKeepsAllTubeVoxels(tc)
            % the Yang-Cheng speed-up mask must never skip a voxel that
            % could respond (lambda2<0 & lambda3<0): compare against
            % eigenvalues computed on every voxel
            rng(2);
            V = TestMfat3D.tube([24 24 24], 2) + 0.2*rand(24,24,24);
            V = V / max(V(:));
            sigma = 1.5;
            geom = mfatCore3D(single(V), sigma, struct('spacing',[1 1 1]));
            [a,b,c,d,e,f] = applyHessian3DAniso(single(V), sigma, [1 1 1]);
            [~, l2, l3] = eig3volume(a,b,c,d,e,f);
            l2(abs(l2) < 1e-4) = 0;  l3(abs(l3) < 1e-4) = 0;
            tube = l2 < 0 & l3 < 0;
            tc.verifyGreaterThan(nnz(tube), 100);
            tc.verifyEqual(geom.lambda2(tube), single(l2(tube)), 'AbsTol', single(1e-5));
            tc.verifyEqual(geom.lambda3(tube), single(l3(tube)), 'AbsTol', single(1e-5));
        end

        function testPrecisionConsistency(tc)
            V = TestMfat3D.tube([24 24 24], 2);
            sig = 1:0.5:2;
            oS = mfatLambda3D(V, sig, 'precision', 'single');
            oD = mfatLambda3D(V, sig, 'precision', 'double');
            tc.verifyLessThan(max(abs(single(oD(:)) - oS(:))), 5e-3);
        end

        function testProbTubeAboveBackground(tc)
            V = TestMfat3D.tube([32 32 32], 2.5);
            p = mfatProb3D(V, 1:0.5:2.5);
            core = V > 0.5;
            tc.verifyGreaterThan(mean(p(core)), mean(p(~imdilate(core, ones(5,5,5)))) + 0.2);
        end

        function testNumericAndCellSigmasAgree(tc)
            % scalar physical sigmas == the same values broadcast per axis
            V = TestMfat3D.tube([24 24 24], 2);
            o1 = mfatLambda3D(V, [1 2]);
            o2 = mfatLambda3D(V, {[1 1 1], [2 2 2]});
            tc.verifyEqual(o1, o2, 'AbsTol', single(1e-6));
        end
    end
end
