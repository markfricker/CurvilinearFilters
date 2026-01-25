function [L1, L2, Vx, Vy] = hessianEigen2D(I, sigma, Precision)
%HESSIANEIGEN2D  Hessian eigenvalues and eigenvectors (2D)
%
%   [L1, L2, Vx, Vy] = hessianEigen2D(I, sigma)
%
%   |L1| <= |L2|
%   (Vx,Vy) is eigenvector corresponding to L2 (normal direction)

if nargin < 3
    Precision = 'single';
end

% Precision
if strcmpi(Precision,'single')
    I = single(I);
else
    I = double(I);
end

% Hessian
[Dxx, Dxy, Dyy] = Hessian2D(I, sigma);

% Scale normalization (Lindeberg)
Dxx = sigma^2 * Dxx;
Dxy = sigma^2 * Dxy;
Dyy = sigma^2 * Dyy;

% Eigen decomposition
[L1, L2, Ix, Iy] = eig2image(Dxx, Dxy, Dyy);

% Return eigenvector of L2
Vx = Ix;
Vy = Iy;
end
