function I = generateJunctionTestImage(N)
% generateJunctionTestImage
% Complex synthetic image for junction recovery testing
%
% Includes:
%   - X, Y, T junctions
%   - overlapping lines
%   - multiple orientations
%   - varying thickness
%   - blobs near junctions
%   - mild noise
%
% Size-safe and deterministic.

if nargin < 1
    N = 256;
end

I = zeros(N,N,'single');

[X,Y] = meshgrid(1:N,1:N);

%% ---- Helper to draw thick lines ----
drawLine = @(x1,y1,x2,y2,w) ...
    abs((y2-y1).*X - (x2-x1).*Y + x2*y1 - y2*x1) ./ ...
    hypot(y2-y1,x2-x1) <= w;

%% ---- Central X junction (thick) ----
I(drawLine(40,40, N-40, N-40, 2)) = 1;
I(drawLine(40,N-40, N-40,40, 2)) = 1;

%% ---- Y junction (mixed thickness) ----
cx = round(0.3*N);
cy = round(0.7*N);
I(drawLine(cx,cy, cx-60,cy-80, 1)) = 1;
I(drawLine(cx,cy, cx+60,cy-80, 2)) = 1;
I(drawLine(cx,cy, cx,cy+80, 3))   = 1;

%% ---- T junction ----
tx = round(0.7*N);
ty = round(0.3*N);
I(drawLine(tx-70,ty, tx+70,ty, 2)) = 1;
I(drawLine(tx,ty, tx,ty+90, 3))    = 1;

%% ---- Parallel overlapping lines ----
for k = -6:6:6
    I(drawLine(30, N/2+k, N-30, N/2+k, 1)) = 1;
end

%% ---- Oblique thin network ----
angles = [15 45 75] * pi/180;
for a = angles
    x1 = round(N/2 - 100*cos(a));
    y1 = round(N/2 - 100*sin(a));
    x2 = round(N/2 + 100*cos(a));
    y2 = round(N/2 + 100*sin(a));
    I(drawLine(x1,y1,x2,y2, 1)) = 1;
end

%% ---- Blobs near junctions ----
blobs = [
    0.30 0.70 6;
    0.50 0.50 4;
    0.70 0.30 5
];

for i = 1:size(blobs,1)
    cx = blobs(i,1)*N;
    cy = blobs(i,2)*N;
    r  = blobs(i,3);
    I = I + exp(-((X-cx).^2 + (Y-cy).^2)/(2*r^2));
end

%% ---- Optical blur + noise ----
I = imgaussfilt(I,1.2);
I = I + 0.03*randn(size(I),'single');

%% ---- Normalise ----
I = I - min(I(:));
I = I / max(I(:));

end
