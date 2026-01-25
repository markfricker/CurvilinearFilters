function I = generateTestImage(N)
% Deterministic Hessian test image (size-safe)

if nargin < 1
    N = 256;
end

I = zeros(N,N,'single');

% --- Relative coordinates ---
x1 = round(0.15*N);
x2 = round(0.75*N);

yThin  = round(0.25*N);
yMed   = round(0.45*N);
yThick = round(0.65*N);

% --- Horizontal lines ---
I(yThin:yThin+1,   x1:x2) = 1;     % thin
I(yMed:yMed+3,     x1:x2) = 1;     % medium
I(yThick:yThick+7, x1:x2) = 1;     % thick

% --- Vertical line ---
I(x1:x2, round(0.8*N):round(0.8*N)+3) = 1;

% --- Diagonal line ---
for k = 1:round(0.5*N)
    r = round(0.25*N) + k;
    c = round(0.25*N) + k;
    if r<=N && c<=N
        I(r,c) = 1;
    end
end

% --- Blobs ---
[X,Y] = meshgrid(1:N,1:N);
blobs = [
    0.30 0.80 0.03;
    0.50 0.80 0.06;
    0.70 0.80 0.10
];

for i = 1:size(blobs,1)
    cx = blobs(i,1)*N;
    cy = blobs(i,2)*N;
    r  = blobs(i,3)*N;
    I = I + exp(-((X-cx).^2 + (Y-cy).^2)/(2*r^2));
end

% --- Optics + noise ---
I = imgaussfilt(I,1);
I = I + 0.02*randn(size(I),'single');

% --- Normalize ---
I = I - min(I(:));
I = I / max(I(:));
end
