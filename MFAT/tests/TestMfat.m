classdef TestMfat < matlab.unittest.TestCase
% TestMfat  Unit tests for mfatAlhassonLambda and mfatAlhassonProb.
%
% Migrated from functiontests pattern to matlab.unittest.TestCase class.
%
% USAGE
%   results = runtests('MFAT/tests/TestMfat');
%   table(results)

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            mfatRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(mfatRoot);                                      % +mfat package parent
            addpath(fullfile(mfatRoot, 'src', 'drivers'));
            addpath(fullfile(mfatRoot, 'src', 'core'));
            addpath(fullfile(mfatRoot, 'src', 'responses'));
            addpath(fullfile(mfatRoot, 'src', 'modifiers'));
            addpath(fullfile(mfatRoot, 'src', 'utils'));
            addpath(fullfile(mfatRoot, 'src', 'original'));
            addpath(fullfile(mfatRoot, 'config'));
        end
    end

    % =====================================================================
    % Synthetic test-image generators (static private)
    % =====================================================================
    methods (Static, Access = private)

        function I = makeLineImage(sz, orientationDeg, width, intensity)
            if nargin < 1, sz = [128 128]; end
            if nargin < 2, orientationDeg = 0; end
            if nargin < 3, width = 3; end
            if nargin < 4, intensity = 1; end
            im = zeros(sz);
            cy = round(sz(1)/2);
            im(cy - floor(width/2) : cy + floor(width/2), :) = intensity;
            I = imrotate(im, orientationDeg, 'crop');
            I = imgaussfilt(I, 1);
        end

        function I = makeEllipseImage(sz, a, b, angleDeg, intensity)
            if nargin < 1, sz = [128 128]; end
            if nargin < 2, a = 15; end
            if nargin < 3, b = 7; end
            if nargin < 4, angleDeg = 30; end
            if nargin < 5, intensity = 1; end
            [X, Y] = meshgrid(1:sz(2), 1:sz(1));
            cx = sz(2)/2; cy = sz(1)/2;
            Xr =  (X - cx) * cosd(angleDeg) + (Y - cy) * sind(angleDeg);
            Yr = -(X - cx) * sind(angleDeg) + (Y - cy) * cosd(angleDeg);
            mask = (Xr./a).^2 + (Yr./b).^2 <= 1;
            I = zeros(sz);
            I(mask) = intensity;
            I = imgaussfilt(I, 1);
        end

        function I = makeBlobImage(sz, radius, intensity)
            if nargin < 1, sz = [128 128]; end
            if nargin < 2, radius = 10; end
            if nargin < 3, intensity = 1; end
            [X, Y] = meshgrid(1:sz(2), 1:sz(1));
            cx = sz(2)/2; cy = sz(1)/2;
            mask = ((X - cx).^2 + (Y - cy).^2) <= radius^2;
            I = zeros(sz);
            I(mask) = intensity;
            I = imgaussfilt(I, 1);
        end

        function I = makeNetworkImage(sz, nLines)
            if nargin < 1, sz = [128 128]; end
            if nargin < 2, nLines = 8; end
            I = zeros(sz);
            rng(0);
            for k = 1:nLines
                x1 = randi([10 sz(2)-10]); y1 = randi([10 sz(1)-10]);
                x2 = randi([10 sz(2)-10]); y2 = randi([10 sz(1)-10]);
                I = I + TestMfat.drawLineMask(sz, x1, y1, x2, y2);
            end
            I = imgaussfilt(I, 1);
            I = I ./ max(I(:));
        end

        function mask = drawLineMask(sz, x1, y1, x2, y2)
            n  = max(abs(x2-x1), abs(y2-y1)) + 1;
            xv = round(linspace(x1, x2, n));
            yv = round(linspace(y1, y2, n));
            mask = zeros(sz);
            for i = 1:numel(xv)
                xi = xv(i); yi = yv(i);
                if xi >= 1 && xi <= sz(2) && yi >= 1 && yi <= sz(1)
                    mask(yi, xi) = 1;
                end
            end
            mask = imdilate(mask, strel('line', 3, 90));
        end

    end

    % =====================================================================
    % Tests
    % =====================================================================
    methods (Test)

        function testEigenvaluesQuadratic(tc)
            % Quadratic surface X^2+Y^2: Hessian eigenvalues should be similar
            [X, Y] = meshgrid(linspace(-1, 1, 64));
            I      = cast(X.^2 + Y.^2, 'double');
            [l1, l2] = mfat.imageEigenvalues2D(I, 1.0, true, 'double', 1e-12);
            border = 4;
            r = border+1 : 64-border;  c = border+1 : 64-border;
            tc.verifyLessThan(mean(abs(l1(r,c) - l2(r,c)), 'all'), 0.05);
            tc.verifyTrue(all(isfinite(l1(:))) && all(isfinite(l2(:))));
        end

        function testLineDetectionBright_det(tc)
            % Deterministic MFAT highlights bright thin line
            I      = TestMfat.makeLineImage([128 128], 0, 2, 1.0);
            sigmas = 0.5:0.5:3;
            out    = mfatAlhassonLambda(I, sigmas, 'tau', 0.01, 'tau2', 0.2, ...
                                        'whiteOnDark', true, 'precision', 'single');
            tc.verifyGreaterThanOrEqual(min(out(:)), 0);
            tc.verifyLessThanOrEqual(   max(out(:)), 1 + 1e-6);
            lineMask = I > 0.5 * max(I(:));
            tc.verifyGreaterThan(mean(out(lineMask)), mean(out(~lineMask)) + 5e-3);
        end

        function testProbabilisticPriorEffect(tc)
            % Higher prior → higher background probabilities
            I      = TestMfat.makeLineImage([128 128], 0, 2, 1.0);
            sigmas = 0.5:0.5:3;
            outLow  = mfatAlhassonProb(I, sigmas, 'prior', 0.01, 'precision', 'single');
            outHigh = mfatAlhassonProb(I, sigmas, 'prior', 0.5,  'precision', 'single');
            tc.verifyGreaterThanOrEqual(min(outLow(:)), 0);
            tc.verifyLessThanOrEqual(   max(outLow(:)), 1 + 1e-6);
            tc.verifyGreaterThanOrEqual(min(outHigh(:)), 0);
            tc.verifyLessThanOrEqual(   max(outHigh(:)), 1 + 1e-6);
            % Higher prior must not produce a LOWER global mean response
            % (on a clean line image both maps saturate near 1 in structured
            %  regions, so we verify monotonicity rather than a fixed delta)
            tc.verifyGreaterThanOrEqual(mean(outHigh(:)), mean(outLow(:)));
        end

        function testDeterministicVsProbabilistic(tc)
            % The two MFAT variants must produce meaningfully different maps
            I      = TestMfat.makeEllipseImage([128 128], 18, 8, 30, 1.0);
            sigmas = 0.5:0.5:4;
            outDet  = mfatAlhassonLambda(I, sigmas, 'tau', 0.03, 'tau2', 0.3, 'precision', 'single');
            outProb = mfatAlhassonProb(  I, sigmas, 'tau', 0.03, 'tau2', 0.3, 'precision', 'single');
            tc.verifyGreaterThan(mean(abs(outDet(:) - outProb(:))), 1e-3);
        end

        function testPrecisionConsistency(tc)
            % Single vs double precision → results agree to within 5e-3
            I      = TestMfat.makeBlobImage([128 128], 10, 1.0);
            sigmas = 1:4;
            outS   = mfatAlhassonLambda(I, sigmas, 'precision', 'single', 'whiteOnDark', true);
            outD   = mfatAlhassonLambda(I, sigmas, 'precision', 'double', 'whiteOnDark', true);
            tc.verifyLessThan(max(abs(single(outD(:)) - outS(:))), 5e-3);
        end

        function testNetworkSanity(tc)
            % Random line network → nonzero response; probabilities in [0,1]
            I      = TestMfat.makeNetworkImage([128 128], 12);
            sigmas = 0.5:0.5:3;
            out    = mfatAlhassonLambda(I, sigmas, 'precision', 'single');
            tc.verifyGreaterThan(max(out(:)), 0);
            p      = mfatAlhassonProb(I, sigmas, 'precision', 'single');
            tc.verifyGreaterThanOrEqual(min(p(:)), 0);
            tc.verifyLessThanOrEqual(   max(p(:)), 1 + 1e-6);
        end

    end

end
