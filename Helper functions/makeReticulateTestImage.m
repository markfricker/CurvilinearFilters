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