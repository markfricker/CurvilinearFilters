classdef TestHessian < matlab.unittest.TestCase
% TestHessian  Unit tests for the Hessian-based curvilinear filters.
%
% Consolidates and migrates the 9 individual functiontests files:
%   test_Hessian2D, test_eig2image, test_hessianEigen2D_scale,
%   test_hessian_direction, test_neuriteness_blob_suppression,
%   test_neuriteness_continuity, test_neuriteness_direction,
%   test_vesselness_basic, test_vesselness_polarity
%
% USAGE
%   results = runtests('Hessian/tests/TestHessian');
%   table(results)
%
% REQUIREMENTS
%   Image Processing Toolbox (imgaussfilt, bwconncomp)

    properties (Constant)
        Tol = 1e-3;
    end

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            hRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(hRoot, 'src'));
            addpath(fullfile(hRoot, 'src', 'engine'));
            addpath(fullfile(hRoot, 'src', 'utils'));
            addpath(fullfile(hRoot, 'src', 'neuriteness'));
            addpath(fullfile(hRoot, 'tests', 'utilities'));
        end
    end

    % =====================================================================
    % Shared test images (static private)
    % =====================================================================
    methods (Static, Access = private)

        function I = lineImage(sz, angleDeg, width)
            % Bright line on dark background
            if nargin < 1, sz = 128; end
            if nargin < 2, angleDeg = 0; end
            if nargin < 3, width = 2; end
            [x, y] = meshgrid(1:sz, 1:sz);
            cx = sz/2; cy = sz/2;
            theta = deg2rad(angleDeg);
            xr = (x - cx) * cos(theta) + (y - cy) * sin(theta);
            I = single(abs(xr) <= width);
            I = imgaussfilt(I, 1);
        end

        function I = blobImage(sz, radius)
            % Centred isotropic Gaussian blob
            if nargin < 1, sz = 64; end
            if nargin < 2, radius = 5; end
            [x, y] = meshgrid(1:sz, 1:sz);
            cx = (sz + 1) / 2; cy = (sz + 1) / 2;
            sig = radius / 2;
            I = single(exp(-((x - cx).^2 + (y - cy).^2) / (2 * sig^2)));
            I = I / max(I(:));
        end

    end

    % =====================================================================
    % applyHessian2D  (migrated from test_Hessian2D.m)
    % =====================================================================
    methods (Test)

        function testHessian2D_constantImage(tc)
            % Constant image → all second derivatives ≈ 0
            I = ones(64);
            [Dxx, Dxy, Dyy] = applyHessian2D(I, 2);
            energy = mean(abs([Dxx(:); Dxy(:); Dyy(:)]));
            tc.verifyLessThan(energy, 0.1);
        end

    end

    % =====================================================================
    % eig2image  (migrated from test_eig2image.m)
    % =====================================================================
    methods (Test)

        function testEig2image_eigenOrdering(tc)
            I = TestHessian.blobImage(64, 5);
            [Dxx, Dxy, Dyy] = applyHessian2D(I, 2);
            [L1, L2, Vx, Vy] = eig2image(Dxx, Dxy, Dyy);
            % |λ1| ≤ |λ2| at every pixel
            tc.verifyTrue(all(abs(L1(:)) <= abs(L2(:)) + eps));
            % Eigenvectors normalised wherever defined
            mag  = hypot(Vx, Vy);
            mask = mag > 0;
            if any(mask(:))
                tc.verifyLessThan(max(abs(mag(mask) - 1)), tc.Tol);
            end
        end

    end

    % =====================================================================
    % hessianEigen2D  (migrated from test_hessianEigen2D_scale.m)
    % =====================================================================
    methods (Test)

        function testHessianEigen2D_scaleNormalisation(tc)
            I = TestHessian.blobImage(64, 4);
            [L1a, ~, ~, ~] = hessianEigen2D(I, 2);
            [L1b, ~, ~, ~] = hessianEigen2D(I, 4);
            ratio = max(abs(L1a(:))) / max(abs(L1b(:)));
            tc.verifyGreaterThan(ratio, 0.2);
            tc.verifyLessThan(ratio, 5);
        end

    end

    % =====================================================================
    % hessian2DFilters direction  (migrated from test_hessian_direction.m)
    % =====================================================================
    methods (Test)

        function testHessian2DFilters_directionAlongLine(tc)
            I = TestHessian.lineImage(128, 0, 2);
            [V, ~, D] = hessian2DFilters(I, ...
                'FilterType', 'vesselness', ...
                'Sigmas', 2, ...
                'Parameters', struct('beta', 0.5, 'c', 15));
            mask = V > 0.5 * max(V(:));
            meanAngle = atan2(mean(sin(D(mask))), mean(cos(D(mask))));
            tc.verifyLessThan(abs(meanAngle), pi/8);
        end

    end

    % =====================================================================
    % neuriteness2D  (migrated from test_neuriteness_*.m)
    % =====================================================================
    methods (Test)

        function testNeuriteness_blobSuppression(tc)
            % Neuriteness should respond to lines far more than blobs
            Iline = TestHessian.lineImage(128, 0, 2);
            Iblob = TestHessian.blobImage(128, 6);
            Nline = neuriteness2D(Iline, 2);
            Nblob = neuriteness2D(Iblob, 2);
            lineSupport = nnz(Nline > 0.5);
            blobSupport = nnz(Nblob > 0.5);
            tc.verifyGreaterThan(lineSupport, 3 * blobSupport);
        end

        function testNeuriteness_continuity(tc)
            % Neuriteness on a tilted line → ≤ 2 connected components
            I       = TestHessian.lineImage(128, 30, 2);
            [N, ~]  = neuriteness2D(I, 2);
            bw      = N > 0.4 * max(N(:));
            cc      = bwconncomp(bw);
            tc.verifyLessThanOrEqual(cc.NumObjects, 2);
        end

        function testNeuriteness_direction(tc)
            % Horizontal line → majority of pixels have direction near 0
            I       = TestHessian.lineImage(128, 0, 2);
            [~, D]  = neuriteness2D(I, 2);
            mask    = abs(D) < pi/4;
            tc.verifyGreaterThan(nnz(mask), 500);
        end

    end

    % =====================================================================
    % hessian2DFilters vesselness  (migrated from test_vesselness_*.m)
    % =====================================================================
    methods (Test)

        function testVesselness_lineDetection(tc)
            I = TestHessian.lineImage(128, 0, 2);
            [V, ~, ~] = hessian2DFilters(I, ...
                'FilterType', 'vesselness', ...
                'Sigmas', 1:4, ...
                'Parameters', struct('beta', 0.5, 'c', 15));
            ridgeVal = mean(V(50:80, 50:80), 'all');
            bgVal    = mean(V(1:20,  1:20),  'all');
            tc.verifyGreaterThan(ridgeVal, 5 * bgVal);
        end

        function testVesselness_polarity(tc)
            % WhiteOnDark=true should greatly outperform WhiteOnDark=false
            % when the line is bright on a dark background
            I = TestHessian.lineImage(128, 0, 2);
            [V1, ~, ~] = hessian2DFilters(I, 'FilterType', 'vesselness', ...
                'Sigmas', 2, 'WhiteOnDark', true);
            [V2, ~, ~] = hessian2DFilters(I, 'FilterType', 'vesselness', ...
                'Sigmas', 2, 'WhiteOnDark', false);
            tc.verifyGreaterThan(max(V1(:)), 5 * max(V2(:)));
        end

    end

end
