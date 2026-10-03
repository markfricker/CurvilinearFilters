classdef TestHessian3D < matlab.unittest.TestCase
% TestHessian3D  Unit tests for the 3D Hessian engine (eig3volume,
% hessianEigen3D, vesselnessResponse3D, platenessResponse3D, hessian3DFilters).
%
% USAGE
%   results = runtests('Hessian/tests/TestHessian3D');
%   table(results)
%
% REQUIREMENTS
%   Image Processing Toolbox (imgaussfilt3)

    properties (Constant)
        Tol = 1e-6;   % eig3volume vs MATLAB eig() -- should be near machine precision
    end

    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            hRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(hRoot, 'src'));
            addpath(fullfile(hRoot, 'src', 'engine'));
        end
    end

    % =========================================================================
    % eig3volume correctness -- against MATLAB's own eig(), per voxel.
    % This is the part of the engine with the least margin for a silent
    % derivation error, so it gets checked against ground truth directly
    % rather than trusted from the algebra alone.
    % =========================================================================
    methods (Test)

        function testEig3volume_matchesBuiltinEig_randomMatrices(tc)
            rng(42);
            n = 200;   % 200 random symmetric 3x3 matrices, as a [n 1 1] "volume"
            a = randn(n,1,1); b = randn(n,1,1); c = randn(n,1,1);
            d = randn(n,1,1); e = randn(n,1,1); f = randn(n,1,1);

            [L1, L2, L3] = eig3volume(a, d, f, b, e, c);

            for i = 1:n
                H = [a(i) d(i) f(i); d(i) b(i) e(i); f(i) e(i) c(i)];
                ref = sort(eig(H));                 % ascending algebraic
                got = sort([L1(i) L2(i) L3(i)]);     % same order for comparison
                tc.verifyEqual(got, ref', 'AbsTol', 1e-8, ...
                    sprintf('mismatch at random matrix %d', i));
            end
        end

        function testEig3volume_magnitudeOrdering(tc)
            % Independent of whether the VALUES match eig() (checked above),
            % confirm the OUTPUT ordering contract (|L1|<=|L2|<=|L3|) holds.
            rng(7);
            n = 500;
            a = randn(n,1,1); b = randn(n,1,1); c = randn(n,1,1);
            d = randn(n,1,1); e = randn(n,1,1); f = randn(n,1,1);
            [L1, L2, L3] = eig3volume(a, d, f, b, e, c);
            tc.verifyTrue(all(abs(L1(:)) <= abs(L2(:)) + 1e-9));
            tc.verifyTrue(all(abs(L2(:)) <= abs(L3(:)) + 1e-9));
        end

        function testEig3volume_diagonalMatrix_exact(tc)
            % Degenerate branch (p1==0, already-diagonal) -- exercised
            % separately since it's a masked special case in the code.
            a = 3; b = -7; c = 1.5;
            z = 0;
            [L1, L2, L3] = eig3volume(a, z, z, b, z, c);
            got = sort(abs([L1 L2 L3]));
            ref = sort(abs([a b c]));
            tc.verifyEqual(got, ref, 'AbsTol', 1e-10);
        end

        function testEig3volume_traceIsPreserved(tc)
            % Sum of eigenvalues must equal the trace, for any matrix --
            % a cheap, independent sanity check that doesn't rely on the
            % sorting/labelling being "correct", just that nothing was
            % dropped or duplicated.
            rng(3);
            n = 100;
            a = randn(n,1,1); b = randn(n,1,1); c = randn(n,1,1);
            d = randn(n,1,1); e = randn(n,1,1); f = randn(n,1,1);
            [L1, L2, L3] = eig3volume(a, d, f, b, e, c);
            tc.verifyEqual(L1+L2+L3, a+b+c, 'AbsTol', 1e-8);
        end
    end

    % =========================================================================
    % Response functions on real synthetic geometry: a tube must score high
    % on vesselness and low on plateness, a sheet the reverse. This is the
    % decisive test -- the whole point of building this engine is that it
    % can tell these apart, which no 2D per-slice method can do at all.
    % =========================================================================
    methods (Test)

        function testVesselness_respondsToTube_notSheet(tc)
            sz = 48;
            tube  = tc.cylinderVolume(sz, 3);
            sheet = tc.slabVolume(sz, 3);

            % one c from both volumes, so the comparison is on a common scale
            pVessel.alpha = 0.5; pVessel.beta = 0.5;
            pVessel.c = max(hessian3DFrangiC(tube, 1:4), hessian3DFrangiC(sheet, 1:4));
            Rtube  = hessian3DFilters(tube,  'FilterType','vesselness', ...
                'Sigmas', 1:4, 'Parameters', pVessel);
            Rsheet = hessian3DFilters(sheet, 'FilterType','vesselness', ...
                'Sigmas', 1:4, 'Parameters', pVessel);

            tubeScore  = max(Rtube(:));
            sheetScore = max(Rsheet(:));
            tc.verifyGreaterThan(tubeScore, 0.01, 'vesselness should respond strongly to a real tube');
            tc.verifyGreaterThan(tubeScore, 3*sheetScore, ...
                'vesselness should respond much more strongly to a tube than a sheet');
        end

        function testPlateness_respondsToSheet_notTube(tc)
            sz = 48;
            tube  = tc.cylinderVolume(sz, 3);
            sheet = tc.slabVolume(sz, 3);

            pPlate.alpha = 0.5; pPlate.beta = 0.5;
            pPlate.c = max(hessian3DFrangiC(tube, 1:4), hessian3DFrangiC(sheet, 1:4));
            Rtube  = hessian3DFilters(tube,  'FilterType','plate', ...
                'Sigmas', 1:4, 'Parameters', pPlate);
            Rsheet = hessian3DFilters(sheet, 'FilterType','plate', ...
                'Sigmas', 1:4, 'Parameters', pPlate);

            % Discrimination is judged at the cores (tube axis, sheet
            % mid-plane). The max over the tube volume sits on its flank
            % (r ~2.5 of 3, sigma 2), where the wall is locally sheet-like
            % (Ra ~0.3) -- real geometry: 0.34 vs 0.86 on the sheet. With
            % the old fixed c = 15 every response was ~1e-4 and ranked by
            % S^2, which hid the weaker-S flank (2026-10-03).
            [X, Y, Z] = ndgrid(1:sz, 1:sz, 1:sz);
            c0 = (sz+1)/2;
            tubeCore  = hypot(X-c0, Y-c0) <= 1;
            sheetCore = abs(Z-c0) <= 1;
            tubeScore  = median(Rtube(tubeCore));
            sheetScore = median(Rsheet(sheetCore));
            tc.verifyGreaterThan(sheetScore, 0.5, 'plateness should respond strongly to a real sheet');
            tc.verifyGreaterThan(sheetScore, 3*tubeScore, ...
                'plateness should respond much more strongly to a sheet than a tube');
            tc.verifyGreaterThan(max(Rsheet(:)), 2*max(Rtube(:)), ...
                'nowhere on the tube should look as sheet-like as the sheet');
        end

        function testPlateness_noisySheet_notGatedByInPlaneSign(tc)
            % A bright sheet curves down only across its thickness (L3);
            % the in-plane L1/L2 are ~0 with a noise-driven sign. The
            % polarity gate used to require L2<0 for 'plate' as well as
            % for vesselness, which zeroed 46% of in-sheet voxels here
            % (sigma=1, 20% noise). Gated on L3 only, ~none should be.
            rng(11);
            n = 64; nz = 25; cz = 13;
            [X, Y, Z] = ndgrid(1:n, 1:n, 1:nz);
            r = hypot(X-(n+1)/2, Y-(n+1)/2);
            I = imgaussfilt3(single(r <= 26 & Z == cz), 1.5);
            I = I + 0.2*max(I(:))*randn(size(I), 'single');
            inSheet = r <= 18 & Z == cz;

            p.alpha = 0.5; p.beta = 0.5; p.c = tc.estimateC(I, 1);
            R = hessian3DFilters(I, 'FilterType','plate', 'Sigmas', 1, 'Parameters', p);
            tc.verifyLessThan(mean(R(inSheet) == 0), 0.02, ...
                'plate response must not be zeroed by the sign of the in-plane eigenvalue L2');

            % Polarity still enforced through L3: a dark sheet on a bright
            % background gives (almost) nothing with WhiteOnDark=true.
            Rdark = hessian3DFilters(max(I(:)) - I, 'FilterType','plate', 'Sigmas', 1, 'Parameters', p);
            tc.verifyLessThan(median(Rdark(inSheet)), 1e-3);
        end

        function testPlateness_suppressesBlob(tc)
            % Descoteaux's blob term |2|L3|-|L2|-|L1||/|L3| is ~0 for a
            % blob (all three eigenvalues equal), so a blob's centre must
            % score far below a sheet's. The previous vesselness-style
            % |L1|/sqrt(|L2*L3|) term gave only ~36x here.
            sz = 48; c0 = (sz+1)/2;
            [X, Y, Z] = ndgrid(1:sz, 1:sz, 1:sz);
            sheet = imgaussfilt3(single(abs(Z-c0) <= 0.5), 1.5);
            blob  = imgaussfilt3(single(sqrt((X-c0).^2+(Y-c0).^2+(Z-c0).^2) <= 2.5), 1.5);
            centre = sqrt((X-c0).^2+(Y-c0).^2+(Z-c0).^2) <= 1;

            p.alpha = 0.5; p.beta = 0.5;
            p.c = tc.estimateC(sheet, 2);
            Rsheet = hessian3DFilters(sheet, 'FilterType','plate', 'Sigmas', 2, 'Parameters', p);
            p.c = tc.estimateC(blob, 2);
            Rblob  = hessian3DFilters(blob,  'FilterType','plate', 'Sigmas', 2, 'Parameters', p);

            tc.verifyGreaterThan(median(Rsheet(centre)), 0.5);
            tc.verifyGreaterThan(median(Rsheet(centre)), 100*median(Rblob(centre)), ...
                'plateness must suppress a blob centre relative to a sheet');
        end

        function testPolarity_darkTubeSuppressedByDefault(tc)
            sz = 32;
            tube = tc.cylinderVolume(sz, 3);   % bright tube on dark background
            darkTube = 1 - tube;               % dark tube on bright background

            p.alpha = 0.5; p.beta = 0.5;       % c from the data (same for both: H(1-f) = -H(f))
            Rbright = hessian3DFilters(tube,     'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p);
            Rdark   = hessian3DFilters(darkTube, 'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p, 'WhiteOnDark', true);
            RdarkOK = hessian3DFilters(darkTube, 'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p, 'WhiteOnDark', false);

            tc.verifyGreaterThan(max(Rbright(:)), 0.01);
            tc.verifyLessThan(max(Rdark(:)), 1e-3, 'WhiteOnDark=true must suppress a dark-on-bright tube');
            tc.verifyGreaterThan(max(Rbright(:)), 10*max(Rdark(:)), ...
                'the bright-tube response must be far stronger than any dark-tube boundary residual');
            tc.verifyGreaterThan(max(RdarkOK(:)), 0.01, 'WhiteOnDark=false must recover the dark-on-bright tube');
        end

        function testIsotropicPath_matchesAnisoUnitSpacing(tc)
            % hessianEigen3D (sigma^2 * applyHessian3D) must equal the
            % anisotropic path at spacing [1 1 1]. Before 2026-10-03
            % applyHessian3D carried the 2D constant 1/(2*pi*sigma^4), so
            % it read sqrt(2*pi)*sigma times too large.
            rng(2);
            V = imgaussfilt3(rand(20, 22, 18), 1);
            for s = [0.7 1.5 3]
                [a1, a2, a3] = hessianEigen3D(V, s, 'double');
                [b1, b2, b3] = hessianEigen3DAniso(V, s, [1 1 1], 'double');
                tol = 1e-12 * max(abs(b3(:)));
                tc.verifyEqual(a1, b1, 'AbsTol', tol);
                tc.verifyEqual(a2, b2, 'AbsTol', tol);
                tc.verifyEqual(a3, b3, 'AbsTol', tol);
            end
        end

        function testScaleNormalisedSheetCurvature_peaksAtMatchedScale(tc)
            % A thin sheet blurred by s0 = 1.5 has a Gaussian profile of
            % width sqrt(s0^2 + sigma^2) after the Hessian's smoothing; its
            % sigma^2-normalised peak curvature sigma^2/(s0^2+sigma^2)^1.5
            % is largest at sigma = sqrt(2)*s0 = 2.12. With the old sigma^3
            % net normalisation it rose monotonically with sigma instead.
            n1 = 81; c1 = 41;
            V = zeros(n1, 9, 9);
            V(c1,:,:) = 1;
            V = imgaussfilt3(V, 1.5);
            sig = 1:0.25:4;
            L3 = zeros(size(sig));
            for i = 1:numel(sig)
                [~, ~, L] = hessianEigen3D(V, sig(i), 'double');
                L3(i) = abs(L(c1,5,5));
            end
            [~, iMax] = max(L3);
            tc.verifyEqual(sig(iMax), sqrt(2)*1.5, 'AbsTol', 0.25);
        end

        function testDefaultC_isIntensityScaleInvariant(tc)
            % With c omitted it comes from the data (hessian3DFrangiC), so
            % scaling the intensity must not change the response. A fixed
            % c (the old default 15) made it depend on the scale.
            tube = tc.cylinderVolume(32, 3);
            p.alpha = 0.5; p.beta = 0.5;
            R1 = hessian3DFilters(tube,       'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p);
            R2 = hessian3DFilters(1000*tube,  'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p);
            tc.verifyGreaterThan(max(R1(:)), 0.1);
            tc.verifyEqual(double(R2), double(R1), 'AbsTol', 1e-4);
        end

        function testScaleSelection_tracksObjectRadius(tc)
            sz = 64;
            thinTube  = tc.cylinderVolume(sz, 2);
            thickTube = tc.cylinderVolume(sz, 6);

            p.alpha = 0.5; p.beta = 0.5;       % c from each volume's own data
            [~, scaleThin]  = hessian3DFilters(thinTube,  'FilterType','vesselness', 'Sigmas', 1:8, 'Parameters', p);
            [~, scaleThick] = hessian3DFilters(thickTube, 'FilterType','vesselness', 'Sigmas', 1:8, 'Parameters', p);

            % Compare the modal selected scale at the tube's own core voxels
            % (where the response is actually large), not the whole volume
            % (mostly background, scale index 0 = "no response anywhere").
            midZ = round(sz/2);
            coreThin  = thinTube(:,:,midZ)  > 0.5;
            coreThick = thickTube(:,:,midZ) > 0.5;
            sThin  = scaleThin(:,:,midZ);
            sThick = scaleThick(:,:,midZ);
            tc.verifyGreaterThan(median(double(sThick(coreThick))), median(double(sThin(coreThin))), ...
                'a thicker tube should select a larger scale index than a thinner one');
        end
    end

    % =========================================================================
    % Synthetic volume helpers
    % =========================================================================
    methods (Access = private, Static)

        function c = estimateC(V, sigma)
            % Frangi c from the data: max Frobenius norm / 2 (as in
            % TestHessian3DAniso.estimateC).
            [L1, L2, L3] = hessianEigen3D(V, sigma, 'single');
            S = sqrt(L1.^2 + L2.^2 + L3.^2);
            c = max(S(:)) / 2;
        end

        function V = cylinderVolume(sz, radius)
            % A bright cylinder (tube) running along Z, centred in X/Y,
            % spanning the full Z extent -- a clean, unambiguous "tube".
            [X, Y, ~] = ndgrid(1:sz, 1:sz, 1:sz);
            cx = (sz+1)/2; cy = (sz+1)/2;
            V = single(hypot(X-cx, Y-cy) <= radius);
            V = imgaussfilt3(V, 1);
        end

        function V = slabVolume(sz, halfThickness)
            % A bright slab (sheet) spanning the full X/Y extent, thin in Z
            % -- a clean, unambiguous "plate".
            [~, ~, Z] = ndgrid(1:sz, 1:sz, 1:sz);
            cz = (sz+1)/2;
            V = single(abs(Z-cz) <= halfThickness);
            V = imgaussfilt3(V, 1);
        end
    end
end
