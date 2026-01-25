%% mfatCompareDemo.m
% Side-by-side demo: deterministic (mfatAlhassonLambda) vs probabilistic (mfatAlhassonProb)
% Displays results for a few synthetic test images. Does NOT save files.

close all;
clearvars;
clc;

% Choose sigmas and common parameters
sigmas = 0.5:0.5:2;
tau = 0.03;
tau2 = 0.3;
precision = 'single';
whiteOnDark = true; % set depending on your synthetic images

% Prepare synthetic images
Iline = makeLineImageDemo([256 256], 0, 2, 1.0);
Iellipse = makeEllipseImageDemo([256 256], 28, 12, 25, 1.0);
Iblob = makeBlobImageDemo([256 256], 18, 1.0);
Inetwork = makeNetworkImageDemo([256 256], 18);

imageSet = {Iline, Iellipse, Iblob, Inetwork};
titles = {'Thin Line', 'Oblong Ellipse', 'Round Blob', 'Simple Network'};

figure('Name','MFAT: Deterministic vs Probabilistic', 'Color', [1 1 1], 'Units','normalized','Position',[0.05 0.05 0.9 0.8]);

for i = 1:numel(imageSet)
    I = imageSet{i};
    % run deterministic
    outDet = mfatAlhassonLambda(I, sigmas, 'tau', tau, 'tau2', tau2, 'D', 0.27, 'whiteOnDark', whiteOnDark, 'precision', precision);
    % run probabilistic
    outProb = mfatAlhassonProb(I, sigmas, 'tau', tau*10, 'tau2', tau2*10, 'whiteOnDark', whiteOnDark, 'prior', 0.1, 'precision', precision, 'finalSmoothSigma', 0.85);

    % Display original, deterministic, probabilistic side-by-side
    idx = (i-1)*3;
    subplot(numel(imageSet), 3, idx + 1);
    imshow(I, []); title([titles{i} ' (orig)']); axis image off;

    subplot(numel(imageSet), 3, idx + 2);
    imshow(outDet, []); title('Deterministic (mfatAlhassonLambda)'); axis image off;
    colormap(gca, gray); colorbar('eastoutside');

    subplot(numel(imageSet), 3, idx + 3);
    imshow(outProb, []); title('Probabilistic (mfatAlhassonProb)'); axis image off;
    colormap(gca, gray); colorbar('eastoutside');
end

% simple helper functions for demo images (local to script)
function I = makeLineImageDemo(sz, orientationDeg, width, intensity)
    if nargin<1, sz=[256 256]; end
    if nargin<2, orientationDeg=0; end
    if nargin<3, width=3; end
    if nargin<4, intensity=1; end
    im = zeros(sz);
    cx = round(sz(2)/2); cy = round(sz(1)/2);
    im(cy-floor(width/2):cy+floor(width/2), :) = intensity;
    I = imrotate(im, orientationDeg, 'crop');
    I = imgaussfilt(I, 1);
end

function I = makeEllipseImageDemo(sz, a, b, angleDeg, intensity)
    if nargin<1, sz=[256 256]; end
    if nargin<2, a=28; end
    if nargin<3, b=12; end
    if nargin<4, angleDeg=25; end
    if nargin<5, intensity=1; end
    [X,Y] = meshgrid(1:sz(2),1:sz(1));
    cx = sz(2)/2; cy = sz(1)/2;
    Xr = ( (X-cx).*cosd(angleDeg) + (Y-cy).*sind(angleDeg) );
    Yr = ( -(X-cx).*sind(angleDeg) + (Y-cy).*cosd(angleDeg) );
    mask = (Xr./a).^2 + (Yr./b).^2 <= 1;
    I = zeros(sz);
    I(mask) = intensity;
    I = imgaussfilt(I, 1);
end

function I = makeBlobImageDemo(sz, radius, intensity)
    if nargin<1, sz=[256 256]; end
    if nargin<2, radius=18; end
    if nargin<3, intensity=1; end
    [X,Y] = meshgrid(1:sz(2),1:sz(1));
    cx = sz(2)/2; cy = sz(1)/2;
    mask = ((X-cx).^2 + (Y-cy).^2) <= radius^2;
    I = zeros(sz); I(mask) = intensity;
    I = imgaussfilt(I, 1);
end

function I = makeNetworkImageDemo(sz, nLines)
    if nargin < 1, sz=[256 256]; end
    if nargin < 2, nLines = 18; end
    I = zeros(sz);
    rng(0);
    for k = 1:nLines
        x1 = randi([10 sz(2)-10]); y1 = randi([10 sz(1)-10]);
        x2 = randi([10 sz(2)-10]); y2 = randi([10 sz(1)-10]);
        mask = drawLineMaskDemo(sz, x1, y1, x2, y2);
        I = I + mask;
    end
    I = imgaussfilt(I, 1);
    I = I ./ max(I(:));
end

function mask = drawLineMaskDemo(sz, x1, y1, x2, y2)
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

% Bring window front
shg;
