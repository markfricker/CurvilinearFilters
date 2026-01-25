function tests = test_phasecong3Optimised
% Unit tests for phasecong3Optimised.m
tests = functiontests(localfunctions);
end

% -------------------------------------------------------------------------
function setupOnce(testCase) %#ok<INUSD>
% Ensure function is on path
assert(exist('phasecong3Optimised','file') == 2, ...
    'phasecong3Optimised.m not found on the MATLAB path.');
end

% -------------------------------------------------------------------------
function testBasicOutputsAndRanges(testCase)

im = makeReticulateTestImage(256, 256, 123); % deterministic synthetic
[M,m,ori,ft,PC,EO,T,pcSum,noiseMask] = phasecong3Optimised( ...
    im, 'precision','double', 'storeEO', true, 'storePC', true);

tol = 1e-8;  % tolerance for numerical roundoff (can be tightened if desired)

% sizes
verifySize(testCase, M, size(im));
verifySize(testCase, m, size(im));
verifySize(testCase, ori, size(im));
verifySize(testCase, ft, size(im));
verifySize(testCase, pcSum, size(im));
verifySize(testCase, noiseMask, size(im));

% finiteness
verifyFalse(testCase, any(~isfinite(M(:))), 'M contains NaN/Inf');
verifyFalse(testCase, any(~isfinite(m(:))), 'm contains NaN/Inf');
verifyFalse(testCase, any(~isfinite(ft(:))), 'featType contains NaN/Inf');
verifyFalse(testCase, any(~isfinite(pcSum(:))), 'pcSum contains NaN/Inf');

% EO shape
verifyTrue(testCase, iscell(EO), 'EO should be a cell when storeEO=true');
verifyEqual(testCase, size(EO,2), 6);        % default norient=6
verifyEqual(testCase, size(EO,1), 4);        % default nscale=4

% PC shape
verifyTrue(testCase, iscell(PC), 'PC should be a cell when storePC=true');
verifyEqual(testCase, numel(PC), 6);

% value ranges / sanity (allow tiny negative due to floating round-off)
verifyGreaterThanOrEqual(testCase, M, -tol);
verifyGreaterThanOrEqual(testCase, m, -tol);

% M should be >= m up to tolerance
verifyGreaterThanOrEqual(testCase, M - m, -tol);

% Orientation range
verifyGreaterThanOrEqual(testCase, ori, 0);
verifyLessThanOrEqual(testCase, ori, 180);

% PC values should be in [0,1] (allow tiny numeric tolerance)
for o = 1:numel(PC)
    verifyFalse(testCase, any(~isfinite(PC{o}(:))), sprintf('PC{%d} contains NaN/Inf', o));
    verifyGreaterThanOrEqual(testCase, PC{o}, -1e-6);
    verifyLessThanOrEqual(testCase, PC{o}, 1+1e-6);
end

% noiseMask logical
verifyTrue(testCase, islogical(noiseMask), 'noiseMask should be logical');

% T scalar
verifyTrue(testCase, isscalar(T), 'T should be a scalar threshold');

end

% -------------------------------------------------------------------------
function testDefaultStoragePolicy(testCase)

im = makeReticulateTestImage(128, 128, 1);

% By default storePC=false (opt-in), EO auto by nargout.
% Request PC/EO in outputs but do not opt in storePC -> PC should be [].
[~,~,~,~,PC,EO] = phasecong3Optimised(im, 'precision','double');

verifyTrue(testCase, isempty(PC), 'PC should be empty by default (storePC=false).');

% EO is only auto-stored if requested as an output.
verifyTrue(testCase, iscell(EO), 'EO should be a cell when requested as output.');

% Now explicitly force no EO
[~,~,~,~,PC2,EO2] = phasecong3Optimised(im, 'precision','double', 'storeEO', false);
verifyTrue(testCase, isempty(EO2), 'EO should be [] when storeEO=false.');
verifyTrue(testCase, isempty(PC2), 'PC should still be [] by default.');

% Explicitly opt-in PC
[~,~,~,~,PC3] = phasecong3Optimised(im, 'precision','double', 'storePC', true);
verifyTrue(testCase, iscell(PC3), 'PC should be a cell when storePC=true.');
end

% -------------------------------------------------------------------------
function testDeterminism(testCase)

im = makeReticulateTestImage(192, 192, 77);

A = phasecong3Optimised(im, 'precision','double', 'storeEO', false, 'storePC', false);
B = phasecong3Optimised(im, 'precision','double', 'storeEO', false, 'storePC', false);

verifyEqual(testCase, A, B, 'AbsTol', 0, ...
    'Outputs should be bitwise identical for deterministic inputs/settings.');
end

% -------------------------------------------------------------------------
function testSingleVsDoubleClose(testCase)

imD = makeReticulateTestImage(192, 192, 9);          % double
imS = single(imD);                                   % single

[M1,m1,or1,ft1,~,~,T1,pcSum1] = phasecong3Optimised( ...
    imD, 'precision','double', 'storeEO', false, 'storePC', false);

[M2,m2,or2,ft2,~,~,T2,pcSum2] = phasecong3Optimised( ...
    imS, 'precision','single', 'storeEO', false, 'storePC', false);

% Expect close, not identical (single precision)
verifyEqual(testCase, M1, double(M2), 'RelTol', 5e-3, 'AbsTol', 5e-4);
verifyEqual(testCase, m1, double(m2), 'RelTol', 5e-3, 'AbsTol', 5e-4);

% Orientation is quantized to integer degrees; should match often, but allow small diffs
verifyEqual(testCase, or1, double(or2), 'AbsTol', 1);

verifyEqual(testCase, ft1, double(ft2), 'RelTol', 5e-3, 'AbsTol', 5e-4);
verifyEqual(testCase, pcSum1, double(pcSum2), 'RelTol', 5e-3, 'AbsTol', 5e-4);

% T may vary slightly
verifyEqual(testCase, T1, double(T2), 'RelTol', 5e-3, 'AbsTol', 5e-4);
end

% -------------------------------------------------------------------------
function testNoiseMaskSanity(testCase)

im = makeReticulateTestImage(128, 128, 101);

% Use fixed threshold so T is known and noiseMask uses Energy_pre <= T.
[~,~,~,~,~,~,T,~,noiseMask] = phasecong3Optimised( ...
    im, 'precision','double', 'noiseMethod', 0.2, 'storeEO', false, 'storePC', false);

verifyTrue(testCase, isscalar(T));
verifyTrue(testCase, islogical(noiseMask));
verifyGreaterThanOrEqual(testCase, nnz(noiseMask), 0);
verifyLessThanOrEqual(testCase, nnz(noiseMask), numel(noiseMask));

% If threshold is very low, very few pixels should be classified as noise
[~,~,~,~,~,~,~,~,noiseMaskLow] = phasecong3Optimised( ...
    im, 'precision','double', 'noiseMethod', 0.0, 'storeEO', false, 'storePC', false);

verifyLessThanOrEqual(testCase, nnz(noiseMaskLow), nnz(noiseMask) + 1); % monotonic-ish sanity
end

% -------------------------------------------------------------------------
function testPersistentCacheLoopUsage(testCase)

im = makeReticulateTestImage(256, 256, 202);

% Call repeatedly; should not error; outputs should stay identical.
ref = phasecong3Optimised(im, 'precision','double', 'storeEO', false, 'storePC', false);

for i = 1:10 %#ok<NASGU>
    out = phasecong3Optimised(im, 'precision','double', 'storeEO', false, 'storePC', false);
    verifyEqual(testCase, out, ref, 'AbsTol', 0);
end

% Optional: clear persistent cache and verify still works and same output
clear phasecong3Optimised
out2 = phasecong3Optimised(im, 'precision','double', 'storeEO', false, 'storePC', false);
verifyEqual(testCase, out2, ref, 'AbsTol', 0);

end

% =========================================================================
% Helper: deterministic synthetic "reticulate-ish" test image (no toolboxes)
% =========================================================================
function im = makeReticulateTestImage(rows, cols, seed)

rng(seed);

% Smooth background via repeated box filtering
im = rand(rows, cols);
im = boxfilt2(im, 7);
im = boxfilt2(im, 11);
im = mat2gray(im)*0.4 + 0.1;

% Add a reticulate line network
nNodes = 40;
pts = [rand(nNodes,1)*(cols-1)+1, rand(nNodes,1)*(rows-1)+1];

% connect each node to a few nearest neighbors
D = squareform(pdist(pts));
D(1:nNodes+1:end) = inf;
[~,idx] = sort(D,2,'ascend');
edges = idx(:,1:3);

canvas = zeros(rows, cols);
for i = 1:nNodes
    for j = edges(i,:)
        canvas = drawThickLine(canvas, pts(i,:), pts(j,:), 1 + randi([0 2]));
    end
end

% soften + add to background
canvas = boxfilt2(canvas, 3);
im = im + 0.8*mat2gray(canvas);

% add noise
im = im + 0.03*randn(rows, cols);

% normalize and cast to double
im = mat2gray(im);
im = double(im);
end

function y = boxfilt2(x, k)
% Simple separable box filter, k odd integer
ker = ones(1,k)/k;
y = conv2(x, ker, 'same');
y = conv2(y, ker', 'same');
end

function img = drawThickLine(img, p0, p1, w)
% Rasterize a thick line by stamping disks along a segment
x0 = p0(1); y0 = p0(2);
x1 = p1(1); y1 = p1(2);
L = hypot(x1-x0, y1-y0);
n = max(2, ceil(L));
xs = linspace(x0, x1, n);
ys = linspace(y0, y1, n);

for t = 1:n
    img = stampDisk(img, xs(t), ys(t), w);
end
end

function img = stampDisk(img, x, y, r)
[rows, cols] = size(img);
cx = round(x); cy = round(y);
x1 = max(1, cx-r); x2 = min(cols, cx+r);
y1 = max(1, cy-r); y2 = min(rows, cy+r);
[X,Y] = meshgrid(x1:x2, y1:y2);
mask = (X-x).^2 + (Y-y).^2 <= r^2;
img(y1:y2, x1:x2) = img(y1:y2, x1:x2) + mask;
end
