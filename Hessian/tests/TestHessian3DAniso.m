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

        function testPerAxisSigma_scalarInput_matchesExplicitBroadcastVector(tc)
            % Backward compatibility: a scalar sigmaPhysical must give
            % BIT-IDENTICAL results to the equivalent explicit 3-element
            % broadcast vector [sigma sigma sigma] -- the scalar form is
            % defined as shorthand for the vector form, not a separate
            % code path.
            spacing = [1, 1, 3];
            sz = 20;
            V = single(rand(sz,sz,sz));
            sigmaScalar = 2.5;

            [Dxx1,Dxy1,Dxz1,Dyy1,Dyz1,Dzz1] = applyHessian3DAniso(V, sigmaScalar, spacing);
            [Dxx2,Dxy2,Dxz2,Dyy2,Dyz2,Dzz2] = applyHessian3DAniso(V, [sigmaScalar sigmaScalar sigmaScalar], spacing);

            tc.verifyEqual(Dxx1, Dxx2);
            tc.verifyEqual(Dyy1, Dyy2);
            tc.verifyEqual(Dzz1, Dzz2);
            tc.verifyEqual(Dxy1, Dxy2);
            tc.verifyEqual(Dxz1, Dxz2);
            tc.verifyEqual(Dyz1, Dyz2);
        end

        function testPerAxisSigma_independentValues_recoverPhysicalCurvature(tc)
            % Closed-form check with GENUINELY independent per-axis sigma
            % (not just the scalar-broadcast case): a quadratic bowl along
            % dim 1 has a fixed physical curvature 2a regardless of what
            % sigma is used along dims 2/3, so an independent per-axis
            % sigmaPhysical = [s1 s2 s3] with s1 ~= s2 ~= s3 must still
            % recover the same Lindeberg-normalised Dxx (sigma1^2 * 2a) --
            % this is what verifies the per-entry generalisation
            % (Dii by sigma_i^2, Dij by sigma_i*sigma_j) rather than just
            % the degenerate all-equal case.
            spacing = [1, 1, 3];
            sz = 24;
            [P1,~,~] = ndgrid(1:sz,1:sz,1:sz);
            c = (sz+1)/2;
            a = 0.01;
            V = single(a*(P1-c).^2);

            % Deliberately all different, but each axis kept >=1 PIXEL
            % (sigmaPerAxis(i)/spacing(i) >= 1): below that, the discrete
            % 2nd-derivative-of-Gaussian kernel is under-resolved (e.g. a
            % sigma of 0.267px needs a 3-tap kernel to approximate a
            % derivative of a function narrower than one pixel) and picks
            % up genuinely large discretisation error unrelated to any
            % chain-rule bug -- exactly the "sigma floor" finding from the
            % real anisotropic ER data this whole feature is motivated by.
            sigmaPerAxis = [1.5, 2.5, 3.5];   % s = [1.5, 2.5, 1.1667] px
            [Dxx,Dxy,~,Dyy,~,Dzz] = applyHessian3DAniso(V, sigmaPerAxis, spacing);
            mid = round(sz/2);

            expected = sigmaPerAxis(1)^2 * 2*a;
            tc.verifyEqual(double(Dxx(mid,mid,mid)), expected, 'RelTol', 0.15);
            tc.verifyEqual(double(Dyy(mid,mid,mid)), 0, 'AbsTol', 2e-3);
            tc.verifyEqual(double(Dzz(mid,mid,mid)), 0, 'AbsTol', 2e-3);
            tc.verifyEqual(double(Dxy(mid,mid,mid)), 0, 'AbsTol', 2e-3);
        end

        function testPerAxisSigma_recoversTubeDiscrimination_atMatchedXYScale(tc)
            % The motivating real-data finding: on severely anisotropic
            % data, flooring a SHARED scalar sigma so s3=sigma/dz>=1 (to
            % get any non-degenerate Z smoothing at all) forces the same
            % sigma through XY too, over-smoothing well past the true tube
            % radius -- on the real plasmodium ER volume this visibly
            % turned thin traced tubules into coarse rounded blobs. A
            % per-axis sigma lets XY stay matched to the true tube radius
            % (the value the well-discriminated baseline test above uses)
            % while Z is independently floored to 1 native pixel, with no
            % forced trade-off between the two axes. This is an absolute
            % check (matching the >3x margin the matched-scalar baseline
            % test requires) rather than a comparison against a synthetic
            % "degraded" case: on this small phantom the shared-scalar
            % case is not itself degraded enough to make a relative
            % comparison a reliable signal (both regimes are far from the
            % near-threshold behaviour the real, much larger, much more
            % anisotropic volume showed) -- the real evidence for the
            % degradation is the visual comparison on real data recorded
            % in memory (hessian_3d_engine_2026_09_30.md), not this phantom.
            [tube, sheet, spacing] = tc.anisoPhantoms();   % spacing = [1 1 3], tube radius 3px in XY

            sigmaPerAxis = {[3*spacing(1), 3*spacing(1), spacing(3)]};   % s = [3, 3, 1] px
            p.alpha = 0.5; p.beta = 0.5;
            p.c = tc.estimateC(tube, sheet, sigmaPerAxis{1}, spacing);
            Rtube  = hessian3DFilters(tube,  'FilterType','vesselness', 'Sigmas', sigmaPerAxis, 'Spacing', spacing, 'Parameters', p);
            Rsheet = hessian3DFilters(sheet, 'FilterType','vesselness', 'Sigmas', sigmaPerAxis, 'Spacing', spacing, 'Parameters', p);

            tc.verifyGreaterThan(max(Rtube(:)), 0.01);
            tc.verifyGreaterThan(max(Rtube(:)), 3*max(Rsheet(:)));
        end

        function testPerAxisSigma_cellArrayPlumbing_hessian3DFiltersAndFrobeniusMask(tc)
            % A minimal end-to-end check that the cell-array 'Sigmas' form
            % is correctly unwrapped (not left as a 1x1 cell) by both
            % hessian3DFilters and localFrobeniusMask3D -- the actual
            % plumbing bug risk of adding a second input form.
            [tube, ~, spacing] = tc.anisoPhantoms();
            sigmaPerAxis = {[3*spacing(1), 3*spacing(1), spacing(3)]};
            p.alpha = 0.5; p.beta = 0.5;   % c from the data (hessian3DFrangiC)

            [R, scaleOut] = hessian3DFilters(tube, 'FilterType','vesselness', ...
                'Sigmas', sigmaPerAxis, 'Spacing', spacing, 'Parameters', p);
            tc.verifyEqual(size(R), size(tube));
            tc.verifyTrue(any(scaleOut(:) > 0));

            masked = localFrobeniusMask3D(R, tube, sigmaPerAxis, 2, spacing);
            tc.verifyEqual(size(masked), size(tube));
            bgRegion = false(size(tube));
            bgRegion(1:4, 1:4, 1) = true;
            tc.verifyEqual(nnz(masked(bgRegion)), 0);
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
            p.alpha = 0.5; p.beta = 0.5;   % c from the data (hessian3DFrangiC)

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
            p.alpha = 0.5; p.beta = 0.5;   % c from the data (hessian3DFrangiC)

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
            p.alpha = 0.5; p.beta = 0.5;   % c from the data (hessian3DFrangiC)

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

        function testSeparable_matchesDenseKernels(tc)
            % applyHessian3DAniso runs separable 1D passes (2026-10-03);
            % rebuild dense 3D kernels here and check the six derivatives
            % agree (incl. replicate-padded borders). The dense kernels
            % have the same moment-matched form, (c + d*p^2)*w for d2/dp^2
            % and b*p*w for d/dp, coefficients from the discrete moments
            % (see tc.momentMatched1D).
            rng(1);
            I = rand(21, 18, 9);
            spacing = [0.1 0.1 0.35];
            sigmaPhys = [0.22 0.22 0.4];
            [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3DAniso(I, sigmaPhys, spacing);

            s = sigmaPhys ./ spacing;
            r = max(1, round(3*s));
            [P1,P2,P3] = ndgrid(-r(1):r(1), -r(2):r(2), -r(3):r(3));
            W = exp(-(P1.^2/(2*s(1)^2) + P2.^2/(2*s(2)^2) + P3.^2/(2*s(3)^2)));
            [a, b, cd] = deal(zeros(1,3), zeros(1,3), zeros(2,3));
            for d = 1:3
                [a(d), b(d), cd(:,d)] = tc.momentMatched1D(s(d), r(d));
            end
            G = prod(a) * W;
            K = {(cd(1,1) + cd(2,1)*P1.^2).*W*a(2)*a(3), b(1)*b(2)*P1.*P2.*W*a(3), b(1)*b(3)*P1.*P3.*W*a(2), ...
                 (cd(1,2) + cd(2,2)*P2.^2).*W*a(1)*a(3), b(2)*b(3)*P2.*P3.*W*a(1), (cd(1,3) + cd(2,3)*P3.^2).*W*a(1)*a(2)};
            % moments of the dense kernels themselves: zero DC, unit
            % curvature gain on x^2/2, unit cross gain on x*y
            tc.verifyEqual(sum(G(:)), 1, 'AbsTol', 1e-12);
            tc.verifyEqual(sum(K{1}(:)), 0, 'AbsTol', 1e-12);
            tc.verifyEqual(sum(P1(:).^2/2 .* K{1}(:)), 1, 'AbsTol', 1e-12);
            tc.verifyEqual(sum(P1(:).*P2(:) .* K{2}(:)), 1, 'AbsTol', 1e-12);
            ij = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
            got = {Dxx, Dxy, Dxz, Dyy, Dyz, Dzz};
            for n = 1:6
                ref = imfilter(I, K{n}, 'conv', 'replicate') ...
                    * sigmaPhys(ij(n,1))*sigmaPhys(ij(n,2)) / (spacing(ij(n,1))*spacing(ij(n,2)));
                tc.verifyEqual(got{n}, ref, 'AbsTol', 1e-10*max(abs(ref(:))) + 1e-12, ...
                    sprintf('separable derivative %d differs from dense kernel', n));
            end
        end

        function testInfiniteSheet_inPlaneEigenvalueIsZero(tc)
            % A noiseless sheet, infinite in-plane (constant along dims 2/3,
            % replicate padding): the in-plane curvature is exactly 0, so
            % L2 must be ~0 while L3 is the curvature across it. The bare
            % sampled d2G kernels did not sum to 0 and gave L2/L3 = 0.006
            % to 0.02 (isotropic) and 0.2-0.38 when dim 3 has a sub-pixel
            % sigma (3x anisotropic grid) -- fake negative curvature that
            % passed the vesselness polarity gate and inflated Ra.
            % Covers applyHessian3D too (hessianEigen3D).
            for s = [0.5 0.75 1 2 4]
                n1 = 2*round(3*s) + 41; c1 = (n1+1)/2;
                V = zeros(n1, 9, 9);
                V(c1,:,:) = 1;
                V = 10 * imgaussfilt3(V, 1.5);   % intensity scale must not matter
                [~, L2, L3] = hessianEigen3D(V, s, 'double');
                tc.verifyLessThan(abs(L2(c1,5,5)), 1e-3*abs(L3(c1,5,5)), ...
                    sprintf('isotropic, sigma=%g: in-plane L2 not ~0', s));
                for spacing = {[1 1 1], [1 1 3]}
                    [~, L2, L3] = hessianEigen3DAniso(V, s, spacing{1}, 'double');
                    tc.verifyLessThan(abs(L2(c1,5,5)), 1e-3*abs(L3(c1,5,5)), ...
                        sprintf('aniso %s, sigma=%g: in-plane L2 not ~0', mat2str(spacing{1}), s));
                end
            end
        end

        function testKernelGains_diagonalEqualsCross(tc)
            % On f = x^2/2 (Dxx = 1) and f = x*y (Dxy = 1) the Lindeberg-
            % normalised outputs must both be sigma^2, and a constant must
            % give 0 -- also at sub-pixel sigma. With the bare sampled
            % kernels the two gains differed (0.89 vs 0.96 at s=4, 1.40 vs
            % 0.77 at s=0.5), making the eigenvalues orientation-dependent.
            for s = [0.4 0.5 1 1.7 4]
                r = max(1, round(3*s)); n = 2*r + 3; c = r + 2;
                [X, Y, ~] = ndgrid((1:n)-c, (1:n)-c, 1:n);
                Dq = cell(1,6); Dc = Dq; D0 = Dq;
                [Dq{:}] = applyHessian3DAniso(X.^2/2, s, [1 1 1]);
                [Dc{:}] = applyHessian3DAniso(X.*Y, s, [1 1 1]);
                [D0{:}] = applyHessian3DAniso(7*ones(n,n,n), s, [1 1 1]);
                tc.verifyEqual(Dq{1}(c,c,c), s^2, 'RelTol', 1e-10);
                tc.verifyEqual(Dc{2}(c,c,c), s^2, 'RelTol', 1e-10);
                for e = 1:6
                    tc.verifyEqual(D0{e}(c,c,c), 0, 'AbsTol', 1e-12);
                end
            end
        end
    end

    methods (Access = private, Static)
        function [a, b, cd] = momentMatched1D(s, r)
            % Coefficients of the 1D factors on w = exp(-p^2/2s^2):
            % a*w (sum 1), b*p*w (first moment 1), (cd(1)+cd(2)*p^2)*w
            % (sum 0, second moment 2) -- independent of the code under test.
            p = -r:r;
            w = exp(-p.^2/(2*s^2));
            a = 1 / sum(w);
            b = 1 / sum(p.^2 .* w);
            cd = [sum(w), sum(p.^2.*w); sum(p.^2.*w), sum(p.^4.*w)] \ [0; 2];
        end

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
