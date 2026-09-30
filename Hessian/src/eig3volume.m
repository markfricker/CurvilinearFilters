function [Lambda1, Lambda2, Lambda3] = eig3volume(Dxx, Dxy, Dxz, Dyy, Dyz, Dzz)
%EIG3VOLUME  Closed-form eigenvalues of a per-voxel real symmetric 3x3 matrix.
%
%   [Lambda1, Lambda2, Lambda3] = eig3volume(Dxx, Dxy, Dxz, Dyy, Dyz, Dzz)
%
% 3D counterpart of eig2image.m. Vectorised, no per-voxel loop and no call
% to MATLAB's eig() (which has no array-valued form) -- uses the standard
% trigonometric closed-form solution for real symmetric 3x3 matrices
% (Smith 1961; see e.g. https://en.wikipedia.org/wiki/Eigenvalue_algorithm
% #3x3_matrices, and the equivalent formulation used by every 3D Frangi/
% Hessian implementation, e.g. Kroon's eig3volume.m for FrangiFilter3D).
%
% Only eigenVALUES are returned -- the 2D engine's eigenVECTOR output
% (orientation) was retired 2026-09-27 because nothing downstream read it
% (every skeleton method ignores enhance-stage orientation); the 3D engine
% mirrors that and never computes eigenvectors at all, avoiding the
% genuinely harder null-space/cross-product step 3D eigenvectors need.
%
% INPUTS
%   Dxx,Dxy,Dxz,Dyy,Dyz,Dzz : per-voxel entries of the symmetric matrix
%                              [Dxx Dxy Dxz; Dxy Dyy Dyz; Dxz Dyz Dzz],
%                              all the same size.
%
% OUTPUTS
%   Lambda1, Lambda2, Lambda3 : eigenvalues sorted so |L1| <= |L2| <= |L3|
%                                (the Frangi/Sato convention used
%                                throughout this codebase), same size as
%                                the inputs.
%
% See also: eig2image, hessianEigen3D, applyHessian3D

classOut = class(Dxx);
a = double(Dxx); b = double(Dyy); c = double(Dzz);
d = double(Dxy); e = double(Dyz); f = double(Dxz);

p1 = d.^2 + f.^2 + e.^2;
q  = (a + b + c) / 3;

% Degenerate (already-diagonal) voxels: off-diagonal terms are all zero,
% so the eigenvalues are simply the diagonal entries. Handled by masking
% rather than branching, so the routine stays fully vectorised.
diagMask = (p1 == 0);

p2 = (a - q).^2 + (b - q).^2 + (c - q).^2 + 2*p1;
p  = sqrt(p2 / 6);
pSafe = p;
pSafe(diagMask) = 1;   % avoid 0/0 on the masked-out voxels below

% B = (1/p) * (A - q*I)
Bxx = (a - q) ./ pSafe;  Byy = (b - q) ./ pSafe;  Bzz = (c - q) ./ pSafe;
Bxy = d ./ pSafe;        Bxz = f ./ pSafe;        Byz = e ./ pSafe;

% r = det(B) / 2
detB = Bxx.*(Byy.*Bzz - Byz.^2) - Bxy.*(Bxy.*Bzz - Byz.*Bxz) + Bxz.*(Bxy.*Byz - Byy.*Bxz);
r = detB / 2;
r = min(max(r, -1), 1);   % numerical safety -- acos domain

phi = acos(r) / 3;

eigA = q + 2*p.*cos(phi);                  % algebraically largest
eigC = q + 2*p.*cos(phi + 2*pi/3);         % algebraically smallest
eigB = 3*q - eigA - eigC;                  % middle (trace = sum of the three)

eigA(diagMask) = a(diagMask);
eigB(diagMask) = b(diagMask);
eigC(diagMask) = c(diagMask);

% Sort each voxel's (eigA,eigB,eigC) triplet ascending by |.|, via a
% 3-element sorting network (3 compare-exchanges) -- fully vectorised,
% no reshape/sub2ind gymnastics needed.
[eigA, eigB] = swapIfGreater(eigA, eigB);
[eigB, eigC] = swapIfGreater(eigB, eigC);
[eigA, eigB] = swapIfGreater(eigA, eigB);

Lambda1 = cast(eigA, classOut);
Lambda2 = cast(eigB, classOut);
Lambda3 = cast(eigC, classOut);
end

% =============================================================================
function [lo, hi] = swapIfGreater(x, y)
%SWAPIFGREATER  Elementwise compare-exchange by |.|: returns (lo,hi) with
% |lo| <= |hi|, swapping x/y wherever |x| > |y|.
swap = abs(x) > abs(y);
lo = x; hi = y;
lo(swap) = y(swap);
hi(swap) = x(swap);
end
