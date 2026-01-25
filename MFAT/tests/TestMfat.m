function tests = TestMfat
% TestMfat  Unit tests for mfatAlhassonLambda and mfatAlhassonProb
tests = functiontests(localfunctions);
end

%% Helper: make simple synthetic images
function I = makeLineImage(sz, orientationDeg, width, intensity)
    if nargin<1, sz = [128 128]; end
    if nargin<2, orientationDeg = 0; end
    if nargin<3, width = 3; end
    if nargin<4, intensity = 1; end
    im = zeros(sz);
    cx = round(sz(2)/2); cy = round(sz(1)/2);
    im(cy-floor(width/2):cy+floor(width/2), :) = intensity;
    I = imrotate(im, orientationDeg, 'crop');
    I = imgaussfilt(I, 1);
end

function I = makeEllipseImage(sz, a, b, angleDeg, intensity)
    if nargin<1, sz = [128 128]; end
    if nargin<2, a = 15; end
    if nargin<3, b = 7; end
    if nargin<4, angleDeg = 30; end
    if nargin<5, intensity = 1; end
    [X,Y] = meshgrid(1:sz(2), 1:sz(1));
    cx = sz(2)/2; cy = sz(1)/2;
    Xr = ( (X-cx).*cosd(angleDeg) + (Y-cy).*sind(angleDeg) );
    Yr = ( -(X-cx).*sind(angleDeg) + (Y-cy).*cosd(angleDeg) );
    mask = (Xr./a).^2 + (Yr./b).^2 <= 1;
    I = zeros(sz);
    I(mask) = intensity;
    I = imgaussfilt(I, 1);
end

function I = makeBlobImage(sz, radius, intensity)
    if nargin<1, sz=[128 128]; end
    if nargin<2, radius=10; end
    if nargin<3, intensity=1; end
    [X,Y] = meshgrid(1:sz(2),1:sz(1));
    cx = sz(2)/2; cy = sz(1)/2;
    mask = ((X-cx).^2 + (Y-cy).^2) <= radius^2;
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
        lineMask = drawLineMask(sz, x1, y1, x2, y2);
        I = I + lineMask;
    end
    I = imgaussfilt(I, 1);
    I = I ./ max(I(:));
end

function mask = drawLineMask(sz, x1, y1, x2, y2)
    % Bresenham-like simple mask using improfile
    n = max(abs(x2-x1), abs(y2-y1)) + 1;
    xv = round(linspace(x1, x2, n));
    yv = round(linspace(y1, y2, n));
    mask = zeros(sz);
    for i = 1:numel(xv)
        xi = xv(i); yi = yv(i);
        if xi>=1 && xi<=sz(2) && yi>=1 && yi<=sz(1)
            mask(yi, xi) = 1;
        end
    end
    mask = imdilate(mask, strel('line', 3, 90));
end

%% Test 1: Eigenvalue helper sanity (quadratic surface)
function testEigenvaluesQuadratic(testCase)
    [X,Y] = meshgrid(linspace(-1,1,64));
    I = X.^2 + Y.^2;
    I = cast(I, 'double');
    [l1, l2] = mfat.imageEigenvalues2D(I, 1.0, true, 'double', 1e-12);
    % eigenvalues should be close to each other (Hessian constant)
    verifyLessThan(testCase, max(abs(l1(:)-l2(:))), 1e-6);
    % all finite
    verifyTrue(testCase, all(isfinite(l1(:))) && all(isfinite(l2(:))));
end

%% Test 2: deterministic MFAT highlights bright thin line (whiteOnDark)
function testLineDetectionBright_det(testCase)
    I = makeLineImage([128 128], 0, 2, 1.0);
    sigmas = 0.5:0.5:3;
    out = mfatAlhassonLambda(I, sigmas, 'tau',0.01,'tau2',0.2,'whiteOnDark', true, 'precision','single');
    verifyGreaterThanOrEqual(testCase, min(out(:)), 0);
    verifyLessThanOrEqual(testCase, max(out(:)), 1 + 1e-6);
    lineMask = I > 0.5*max(I(:));
    meanLine = mean(out(lineMask));
    meanBg = mean(out(~lineMask));
    verifyGreaterThan(testCase, meanLine, meanBg + 5e-3);
end

%% Test 3: probabilistic MFAT in [0,1], prior affects mean
function testProbabilisticPriorEffect(testCase)
    I = makeLineImage([128 128], 0, 2, 1.0);
    sigmas = 0.5:0.5:3;
    outLowPrior = mfatAlhassonProb(I, sigmas, 'prior', 0.01, 'precision', 'single');
    outHighPrior = mfatAlhassonProb(I, sigmas, 'prior', 0.5, 'precision', 'single');
    verifyGreaterThanOrEqual(testCase, min(outLowPrior(:)), 0);
    verifyLessThanOrEqual(testCase, max(outLowPrior(:)), 1 + 1e-6);
    verifyGreaterThan(testCase, mean(outHighPrior(:)), mean(outLowPrior(:)) + 1e-3);
end

%% Test 4: deterministic vs probabilistic outputs differ (sanity)
function testDeterministicVsProbabilistic(testCase)
    I = makeEllipseImage([128 128], 18, 8, 30, 1.0);
    sigmas = 0.5:0.5:4;
    outDet = mfatAlhassonLambda(I, sigmas, 'tau',0.03,'tau2',0.3,'precision','single');
    outProb = mfatAlhassonProb(I, sigmas, 'tau',0.03,'tau2',0.3,'precision','single');
    % maps should not be identical; check L1 difference
    diffL1 = mean(abs(outDet(:) - outProb(:)));
    verifyGreaterThan(testCase, diffL1, 1e-3);
end

%% Test 5: precision consistency single vs double
function testPrecisionConsistency(testCase)
    I = makeBlobImage([128 128], 10, 1.0);
    sigmas = 1:1:4;
    outS = mfatAlhassonLambda(I, sigmas, 'precision', 'single', 'whiteOnDark', true);
    outD = mfatAlhassonLambda(I, sigmas, 'precision', 'double', 'whiteOnDark', true);
    diff = single(outD) - single(outS);
    verifyLessThan(testCase, max(abs(diff(:))), 5e-3);
end

%% Test 6: network detection sanity (no crash + nonzero response)
function testNetworkSanity(testCase)
    I = makeNetworkImage([128 128], 12);
    sigmas = 0.5:0.5:3;
    out = mfatAlhassonLambda(I, sigmas, 'precision', 'single');
    verifyGreaterThan(testCase, max(out(:)), 0);
    % probabilistic returns probabilities in [0,1]
    p = mfatAlhassonProb(I, sigmas, 'precision', 'single');
    verifyGreaterThanOrEqual(testCase, min(p(:)), 0);
    verifyLessThanOrEqual(testCase, max(p(:)), 1 + 1e-6);
end
