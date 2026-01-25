function [Lambda1, Lambda2, Ix, Iy] = eig2image(Dxx, Dxy, Dyy)
%EIG2IMAGE Fast eigenvalues and eigenvectors of 2x2 Hessian matrix
%
%   [Lambda1, Lambda2, Ix, Iy] = eig2image(Dxx, Dxy, Dyy)
%
%   Eigenvalues are sorted such that:
%       |Lambda1| <= |Lambda2|
%
%   (Ix, Iy) is the eigenvector corresponding to Lambda2
%   (minor curvature / vessel direction)

% --- eigenvalues ---

tmp = hypot(Dxx - Dyy, 2*Dxy);

mu1 = 0.5 * (Dxx + Dyy + tmp);
mu2 = 0.5 * (Dxx + Dyy - tmp);

% --- sort by absolute value ---
swap = abs(mu1) > abs(mu2);

Lambda1 = mu1;
Lambda2 = mu2;

Lambda1(swap) = mu2(swap);
Lambda2(swap) = mu1(swap);

% --- eigenvector for Lambda2 ---
% Solve (H - Lambda2*I)v = 0
Ix = 2 * Dxy;
Iy = Dyy - Dxx + tmp;

% Handle swapped cases
Ix(swap) = 2 * Dxy(swap);
Iy(swap) = Dyy(swap) - Dxx(swap) - tmp(swap);

% --- normalise eigenvectors ---
mag = hypot(Ix, Iy);
mag(mag == 0) = 1;

Ix = Ix ./ mag;
Iy = Iy ./ mag;
end

