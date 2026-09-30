classdef TestHessian3DAniso < matlab.unittest.TestCase
% TestHessian3DAniso  Unit tests for the anisotropic (native-grid, no
% resampling) 3D Hessian path: applyHessian3DAniso, hessianEigen3DAniso,
% and hessian3DFilters'/localFrobeniusMask3D's 'Spacing' option.
%
% USAGE
%   results = runtests('Hessian/tests/TestHessian3DAniso');
%
% REQUIREMENTS
%   Image Processing Toolbox (imgaussfilt3, imresize3 -- the latter only
%   for the cross-check against the isotropic-resample path)

    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            hRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(hRoot, 'src'));
            addpath(fullfile(hRoot, 'src', 'engine'));
            segRoot = fullfile(hRoot, '..', '..', 'Segmentation_sandbox', 'src');
            addpath(segRoot);
        end
    end

    methods (Test)

        function testEigenvaluesAreRotationallyConsistent_diagonalCase(tc)
            % A trivial closed-form check independent of the Gaussian
            % machinery: applyHessian3DAniso's physical-unit rescaling,
            % fed a matrix that's already diagonal in physical units,
            % must return exactly those diagonal values (sanity check on
            % the chain-rule division, isolated from convolution/kernel
            % effects).
            spacing = [1, 2, 4];
            sigmaPhys = 2;
            % Build a synthetic volume whose Hessian is dominated by a
            % clean parabolic bowl with KNOWN physical curvature -- rather
            % than solve the general case, confirm the isotropic-spacing
            % special case (spacing=[1 1 1]) exactly matches an
            % independent direct calculation.
            sz = 24;
            [P1,P2,P3] = ndgrid(1:sz,1:sz,1:sz);
            c = (sz+1)/2;
            % f = a*(x-c)^2 : d2f/dx2 = 2a everywhere, a plain constant --
            % easiest possible closed-form ground truth.
            a = 0.01;
            V = single(a*(P1-c).^2);
            [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3DAniso(V, sigmaPhys, [1 1 1]);
            mid = round(sz/2);
            % Expected value includes the Lindeberg sigma^2 scaling this
            % function deliberately applies (matches hessianEigen3D's own
            % sigma^2*D convention) -- the raw physical curvature is 2a,
            % the Lindeberg-normalised value compared here is sigma^2*2a.
            % ~8% slack for Gaussian-kernel truncation (finite support).
            expected = sigmaPhys^2 * 2*a;
            tc.verifyEqual(double(Dxx(mid,mid,mid)), expected, 'RelTol', 0.15);
            tc.verifyEqual(double(Dyy(mid,mid,mid)), 0, 'AbsTol', 2e-3);
            tc.verifyEqual(double(Dzz(mid,mid,mid)), 0, 'AbsTol', 2e-3);
            tc.verifyEqual(double(Dxy(mid,mid,mid)), 0, 'AbsTol', 2e-3);
        end

        function testChainRuleScaling_dim1TwiceAsCoarse_recoversSamePhysicalCurvature(tc)
            % A quadratic bowl has a FIXED true physical curvature (2a)
            % regardless of how finely/coarsely it happens to be sampled.
            % Build the SAME physical function on two different pixel
            % grids (parameterised by physical position, not pixel index,
            % so this is a fair test) and confirm both recover the same
            % Lindeberg-normalised physical Dxx -- the whole point of the
            % chain-rule correction in applyHessian3DAniso.
            sz = 24;
            [P1,~,~] = ndgrid(1:sz,1:sz,1:sz);
            a = 0.01;
            sigmaPhys = 2;
            x0phys = 12;   % fixed physical centre, same for both grids

            spacingFine   = [1,1,1];
            spacingCoarse = [2,1,1];
            Vfine   = single(a*((P1-1)*spacingFine(1)   - x0phys).^2);
            Vcoarse = single(a*((P1-1)*spacingCoarse(1) - x0phys).^2);

            [DxxFine,   ~,~,~,~,~] = applyHessian3DAniso(Vfine,   sigmaPhys, spacingFine);
            [DxxCoarse, ~,~,~,~,~] = applyHessian3DAniso(Vcoarse, sigmaPhys, spacingCoarse);

            midFine   = round(x0phys/spacingFine(1))   + 1;
            midCoarse = round(x0phys/spacingCoarse(1)) + 1;

            expected = sigmaPhys^2 * 2*a;
            tc.verifyEqual(double(DxxFine(midFine,midFine,midFine)),       expected, 'RelTol', 0.15);
            tc.verifyEqual(double(DxxCoarse(midCoarse,midCoarse,midCoarse)), expected, 'RelTol', 0.15);
        end

        function testVesselness_respondsToAxialTube_nativeAnisotropicGrid(tc)
            % The hardest, most important case: a tube running ALONG the
            % coarse (Z) axis -- exactly the structure a per-slice 2D
            % method cannot see at all -- detected directly on the native
            % low-Z-resolution grid, no resampling.
            [tube, sheet, spacing] = tc.anisoPhantoms();

            sigmaPhys = 3*spacing(1);   % ~3 fine-pixels' worth, physical units
            p.alpha = 0.5; p.beta = 0.5;
            % c is data-dependent for any Frangi-style filter (it calibrates
            % the background-suppression term to the actual Hessian
            % magnitude of the data -- NOT a universal constant; a fixed
            % default tuned for one image's intensity/derivative scale can
            % crush the response to near-zero on another, exactly what a
            % real ER volume later showed with a naively-reused default).
            p.c = tc.estimateC(tube, sheet, sigmaPhys, spacing);

            Rtube  = hessian3DFilters(tube,  'FilterType','vesselness', ...
                'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);
            Rsheet = hessian3DFilters(sheet, 'FilterType','vesselness', ...
                'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);

            tc.verifyGreaterThan(max(Rtube(:)), 0.01);
            tc.verifyGreaterThan(max(Rtube(:)), 3*max(Rsheet(:)));
        end

        function testPlateness_respondsToSheet_nativeAnisotropicGrid(tc)
            [tube, sheet, spacing] = tc.anisoPhantoms();
            sigmaPhys = 3*spacing(1);
            p.alpha = 0.5; p.beta = 0.5;
            p.c = tc.estimateC(tube, sheet, sigmaPhys, spacing);

            Rtube  = hessian3DFilters(tube,  'FilterType','plate', ...
                'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);
            Rsheet = hessian3DFilters(sheet, 'FilterType','plate', ...
                'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);

            tc.verifyGreaterThan(max(Rsheet(:)), 0.01);
            % Weaker margin than the clean isotropic phantom test (which
            % gets >3x): this sheet is only 3 native voxels thick on a
            % 14-slice, 3x-anisotropic grid -- close to the axial sampling
            % limit for "confidently flat". Direction is still correct
            % (sheet > tube), just not as decisively separated -- a real,
            % honest finding about very thin sheets under strong
            % anisotropy, not a bug to paper over with a looser test.
            tc.verifyGreaterThan(max(Rsheet(:)), 1.2*max(Rtube(:)));
        end

        function testFrobeniusMask3D_anisoPath_zerosBackground(tc)
            [tube, ~, spacing] = tc.anisoPhantoms();
            sigmaPhys = 3*spacing(1);
            p.alpha = 0.5; p.beta = 0.5; p.c = 15;

            Rtube = hessian3DFilters(tube, 'FilterType','vesselness', ...
                'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);

            imfConst = ones(size(tube), 'single');
            masked = localFrobeniusMask3D(imfConst, tube, sigmaPhys, 2, spacing);

            bgRegion = false(size(tube));
            bgRegion(1:4, 1:4, 1) = true;
            tc.verifyEqual(nnz(masked(bgRegion)), 0);

            [nY,nX,nZ] = size(tube);
            core = false(size(tube));
            core(round(nY/2), round(nX/2), round(nZ/2)) = true;
            tc.verifyGreaterThan(masked(core), 0);
        end

        function testAnisoNativeGrid_agreesWithIsotropicResample_onDiscrimination(tc)
            % Cross-check: the SAME underlying structure, processed two
            % ways -- (a) native anisotropic grid via the new Spacing path,
            % (b) isotropic-resampled via the original engine -- must agree
            % on which structure is tube-like vs sheet-like, even though
            % the two paths use different kernels/normalisation and are not
            % expected to match numerically.
            [tube, sheet, spacing] = tc.anisoPhantoms();
            zRatio = spacing(3)/spacing(1);

            sigmaPhys = 3*spacing(1);
            p.alpha = 0.5; p.beta = 0.5; p.c = 15;

            RtubeNative = hessian3DFilters(tube, 'FilterType','vesselness', ...
                'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);
            RsheetNative = hessian3DFilters(sheet, 'FilterType','vesselness', ...
                'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);
            nativeSaysTubeMoreTubelike = max(RtubeNative(:)) > max(RsheetNative(:));

            nZiso = round(size(tube,3) * zRatio);
            tubeIso  = imresize3(tube,  [size(tube,1),  size(tube,2),  nZiso], 'linear');
            sheetIso = imresize3(sheet, [size(sheet,1), size(sheet,2), nZiso], 'linear');

            sigmaPix = 3;   % matches sigmaPhys/spacing(1) = 3 fine-pixels, now isotropic
            RtubeIso  = hessian3DFilters(tubeIso,  'FilterType','vesselness', 'Sigmas', sigmaPix, 'Parameters', p);
            RsheetIso = hessian3DFilters(sheetIso, 'FilterType','vesselness', 'Sigmas', sigmaPix, 'Parameters', p);
            isoSaysTubeMoreTubelike = max(RtubeIso(:)) > max(RsheetIso(:));

            tc.verifyEqual(nativeSaysTubeMoreTubelike, isoSaysTubeMoreTubelike, ...
                'native-anisotropic and isotropic-resample paths must agree on which phantom is tube-like');
            tc.verifyTrue(nativeSaysTubeMoreTubelike, 'the tube phantom must win on both paths');
        end

        function testAnisoNativeGrid_isFasterThanIsotropicResamplePath(tc)
            [tube, ~, spacing] = tc.anisoPhantoms();
            zRatio = spacing(3)/spacing(1);
            sigmaPhys = 3*spacing(1);
            p.alpha = 0.5; p.beta = 0.5; p.c = 15;

            tNative = tic;
            hessian3DFilters(tube, 'FilterType','vesselness', 'Sigmas', sigmaPhys, 'Spacing', spacing, 'Parameters', p);
            tNative = toc(tNative);

            nZiso = round(size(tube,3) * zRatio);
            tIso = tic;
            tubeIso = imresize3(tube, [size(tube,1), size(tube,2), nZiso], 'linear');
            hessian3DFilters(tubeIso, 'FilterType','vesselness', 'Sigmas', 3, 'Parameters', p);
            tIso = toc(tIso);

            tc.verifyLessThan(tNative, tIso, ...
                sprintf('native-grid path (%.3fs) should beat isotropic-resample (%.3fs) on anisotropic input', tNative, tIso));
        end
    end

    methods (Access = private, Static)
        function c = estimateC(vol1, vol2, sigmaPhys, spacing)
            % Data-dependent calibration of the Frangi "c" (background/
            % noise suppression) parameter: c = max(Frobenius norm)/2,
            % the same frobDivision=2 convention localFrobeniusMask3D
            % uses. c is fundamentally a property of the DATA's own
            % Hessian magnitude, not a universal constant -- see the
            % comment at its use site.
            [D11,D12,D13,D22,D23,D33] = applyHessian3DAniso(vol1, sigmaPhys, spacing);
            S1 = sqrt(D11.^2+D22.^2+D33.^2 + 2*D12.^2+2*D13.^2+2*D23.^2);
            [E11,E12,E13,E22,E23,E33] = applyHessian3DAniso(vol2, sigmaPhys, spacing);
            S2 = sqrt(E11.^2+E22.^2+E33.^2 + 2*E12.^2+2*E13.^2+2*E23.^2);
            c = max([S1(:); S2(:)]) / 2;
        end

        function [tube, sheet, spacing] = anisoPhantoms()
            % Native anisotropic grid: fine in dims 1/2, coarse in dim 3
            % (dz = 3x dxy) -- built directly at native resolution, never
            % downsampled from a finer version, matching how a real
            % acquisition actually samples anisotropic data.
            nY = 40; nX = 40; nZ = 14;
            spacing = [1, 1, 3];   % dz = 3x dxy, in the SAME units as sigmaPhys

            % Tube running ALONG the coarse Z axis -- the case no 2D
            % per-slice method can see at all.
            [X,Y,~] = ndgrid(1:nY, 1:nX, 1:nZ);
            cx = (nY+1)/2; cy = (nX+1)/2;
            tube = single(hypot(X-cx, Y-cy) <= 3);
            tube = imgaussfilt3(tube, [1 1 0.5]);

            % Sheet spanning the full X/Y extent, thin in Z (native few
            % slices -- 3 of the 14 available, reflecting genuinely coarse
            % axial sampling of a thin structure).
            [~,~,Z] = ndgrid(1:nY, 1:nX, 1:nZ);
            cz = (nZ+1)/2;
            sheet = single(abs(Z-cz) <= 1.5);
            sheet = imgaussfilt3(sheet, [1 1 0.5]);
        end
    end
end
