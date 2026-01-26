function demo_generate_complex_phantom()
%DEMO_GENERATE_COMPLEX_PHANTOM
% Generates a synthetic benchmark phantom and saves it as a PNG.
%
% Structures:
%  - thin and thick lines
%  - multiple-angle intersections
%  - blobs
%  - plate-like regions

sz = 512;
I = zeros(sz);

[x,y] = meshgrid(1:sz,1:sz);
cx = sz/2; cy = sz/2;

% -------------------------------------------------------------------------
% Thin lines
% -------------------------------------------------------------------------
theta = deg2rad(20);
xr = (x-cx)*cos(theta) + (y-cy)*sin(theta);
I(abs(xr) < 2) = 1;

theta = deg2rad(-40);
xr = (x-cx)*cos(theta) + (y-cy)*sin(theta);
I(abs(xr) < 1.5) = 1;

% -------------------------------------------------------------------------
% Thick line
% -------------------------------------------------------------------------
theta = deg2rad(75);
xr = (x-cx)*cos(theta) + (y-cy)*sin(theta);
I(abs(xr) < 6) = 1;

% -------------------------------------------------------------------------
% Junction (cross)
% -------------------------------------------------------------------------
I(abs(x-cx) < 2 | abs(y-cy) < 2) = 1;

% -------------------------------------------------------------------------
% Blobs
% -------------------------------------------------------------------------
blobs = [
    100 400 12
    400 100 18
    420 420 10
];
for k = 1:size(blobs,1)
    bx = blobs(k,1);
    by = blobs(k,2);
    r  = blobs(k,3);
    I = I + exp(-((x-bx).^2 + (y-by).^2)/(2*r^2));
end

% -------------------------------------------------------------------------
% Plate-like region
% -------------------------------------------------------------------------
I(300:380, 200:420) = 1;

% -------------------------------------------------------------------------
% Smooth and normalize
% -------------------------------------------------------------------------
I = imgaussfilt(I,1);
I = I - min(I(:));
I = I / max(I(:));

% -------------------------------------------------------------------------
% Save phantom
% -------------------------------------------------------------------------
imwrite(I, '../manual/figures/complex_phantom.png');

% Optional visualization
figure; imagesc(I); axis image off; colormap gray;
title('Complex synthetic phantom');
end
