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

            pVessel.alpha = 0.5; pVessel.beta = 0.5; pVessel.c = 15;
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

            pPlate.alpha = 0.5; pPlate.beta = 0.5; pPlate.c = 15;
            Rtube  = hessian3DFilters(tube,  'FilterType','plate', ...
                'Sigmas', 1:4, 'Parameters', pPlate);
            Rsheet = hessian3DFilters(sheet, 'FilterType','plate', ...
                'Sigmas', 1:4, 'Parameters', pPlate);

            tubeScore  = max(Rtube(:));
            sheetScore = max(Rsheet(:));
            tc.verifyGreaterThan(sheetScore, 0.01, 'plateness should respond strongly to a real sheet');
            tc.verifyGreaterThan(sheetScore, 3*tubeScore, ...
                'plateness should respond much more strongly to a sheet than a tube');
        end

        function testPolarity_darkTubeSuppressedByDefault(tc)
            sz = 32;
            tube = tc.cylinderVolume(sz, 3);   % bright tube on dark background
            darkTube = 1 - tube;               % dark tube on bright background

            p.alpha = 0.5; p.beta = 0.5; p.c = 15;
            Rbright = hessian3DFilters(tube,     'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p);
            Rdark   = hessian3DFilters(darkTube, 'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p, 'WhiteOnDark', true);
            RdarkOK = hessian3DFilters(darkTube, 'FilterType','vesselness', 'Sigmas', 1:3, 'Parameters', p, 'WhiteOnDark', false);

            tc.verifyGreaterThan(max(Rbright(:)), 0.01);
            tc.verifyLessThan(max(Rdark(:)), 1e-3, 'WhiteOnDark=true must suppress a dark-on-bright tube');
            tc.verifyGreaterThan(max(Rbright(:)), 10*max(Rdark(:)), ...
                'the bright-tube response must be far stronger than any dark-tube boundary residual');
            tc.verifyGreaterThan(max(RdarkOK(:)), 0.01, 'WhiteOnDark=false must recover the dark-on-bright tube');
        end

        function testScaleSelection_tracksObjectRadius(tc)
            sz = 64;
            thinTube  = tc.cylinderVolume(sz, 2);
            thickTube = tc.cylinderVolume(sz, 6);

            p.alpha = 0.5; p.beta = 0.5; p.c = 15;
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
