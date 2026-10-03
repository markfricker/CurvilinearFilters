function geom = mfatCore3D(V, sigma, opts)
% =========================================================================
% MFAT – Multiscale Fractional Anisotropy Tensor Framework
%
% Author:
%   H. Alhasson, M. Alharbi, B. Obara
%
% Refactored framework & extensions:
%   MD Fricker, Jan 2026 (2D); 3D port Oct 2026
%
% Citation:
%   H. Alhasson, M. Alharbi, B. Obara,
%   "2D and 3D Vascular Structures Enhancement via
%    Multiscale Fractional Anisotropy Tensor",
%   ECCV Workshops (BioImage Computing), 2018.
%
% License:
%   Academic / research use. Please cite the above work.
% =========================================================================

% MFATCORE3D  MFAT geometry core for a volume (single scale)
%
% OVERVIEW
%   3D counterpart of mfatCore2D. Computes, at one scale:
%     - the two LARGEST-magnitude Hessian eigenvalues lambda2, lambda3
%       (|lambda1| <= |lambda2| <= |lambda3|; lambda1, the along-tube
%       eigenvalue, is not used -- as in the published 3D filter),
%     - the three-factor regularisation of lambda3 (tau and tau2 clips),
%     - the fractional-anisotropy-like measure of (|l2|, |l3tau|, |l4|).
%   Purely geometric: no thresholding, no aggregation.
%
% DIFFERENCES FROM THE PUBLISHED 3D CODE (FractionalIstropicTensor3D.m,
% github.com/Haifafh/MFAT) -- deliberate:
%   1. Smoothing is genuinely 3D. The published imgaussian() smooths only
%      along dims 1 and 2, so the Z derivatives there are raw finite
%      differences of unsmoothed data.
%   2. Anisotropic voxels are handled on the native grid (per-axis pixel
%      sigma + physical-unit, Lindeberg-normalised Hessian, via
%      applyHessian3DAniso). The published code takes one scalar spacing
%      and applies it to X and Y only.
%   3. Gaussian-derivative kernels (separable) instead of
%      smooth-then-central-difference-twice, which widens and biases the
%      stencil.
%   4. The tau reference min(lambda3) is taken over tube-candidate voxels
%      (lambda2<0 & lambda3<0, the only voxels that can respond), so the
%      result does not depend on which voxels the Yang-Cheng speed-up mask
%      happened to skip.
%   5. Pure MATLAB eigenvalues (eig3volume, vectorised closed form) -- the
%      published code needs a compiled MEX (Windows/macOS binaries only).
%   6. Single precision throughout (the published code allocates the
%      eigenvalue volumes in double).
%
% INPUTS
%   V     - 3D volume (already normalised).
%   sigma - Gaussian scale in PHYSICAL units: scalar, or [s1 s2 s3] per
%           axis (see applyHessian3DAniso).
%   opts  - Struct with fields .tau, .tau2, .whiteOnDark, .precision,
%           .spacing ([1 1 1] = plain pixel units).
%
% OUTPUT
%   geom - Struct: .lambda2, .lambda3 (raw), .lambda3t (tau clip),
%          .lambda4 (tau2 clip), .fa
%
% See also: mfatCore2D, mfatResponseLambda3D, applyHessian3DAniso, eig3volume

% ---- defaults ----
if ~isfield(opts,'whiteOnDark'), opts.whiteOnDark = true;     end
if ~isfield(opts,'tau'),         opts.tau = 0.03;             end
if ~isfield(opts,'tau2'),        opts.tau2 = 0.3;             end
if ~isfield(opts,'precision'),   opts.precision = 'single';   end
if ~isfield(opts,'spacing'),     opts.spacing = [1 1 1];      end

% ---- precision ----
if strcmpi(opts.precision,'single')
    rc = 'single';
    eps0 = eps('single');
    eigTol = single(1e-4);
else
    rc = 'double';
    eps0 = eps('double');
    eigTol = 1e-8;
end

V = cast(V, rc);

% ---- (1) Hessian (physical units, Lindeberg-normalised) ----
[Hxx,Hxy,Hxz,Hyy,Hyz,Hzz] = applyHessian3DAniso(V, sigma, opts.spacing);

% sign convention: bright structures -> negative eigenvalues
if ~opts.whiteOnDark
    Hxx = -Hxx; Hxy = -Hxy; Hxz = -Hxz;
    Hyy = -Hyy; Hyz = -Hyz; Hzz = -Hzz;
end

% ---- (2) Yang & Cheng (2014) candidate mask ----
% Characteristic polynomial l^3 + B1 l^2 + B2 l + B3; the excluded voxels
% cannot have two negative eigenvalues, so they cannot respond.
B1 = -(Hxx + Hyy + Hzz);
B2 = Hxx.*Hyy + Hxx.*Hzz + Hyy.*Hzz - Hxy.^2 - Hxz.^2 - Hyz.^2;
B3 = Hxx.*Hyz.^2 + Hxy.^2.*Hzz + Hxz.^2.*Hyy - Hxx.*Hyy.*Hzz - 2*Hxy.*Hyz.*Hxz;
T = B1 > 0 & ~(B2 <= 0 & B3 == 0) & ~(B2 > 0 & B1.*B2 < B3);
clear B1 B2 B3
idx = find(T);
clear T

lambda2 = zeros(size(V), rc);
lambda3 = zeros(size(V), rc);
if ~isempty(idx)
    [~, l2, l3] = eig3volume(Hxx(idx), Hxy(idx), Hxz(idx), Hyy(idx), Hyz(idx), Hzz(idx));
    lambda2(idx) = cast(l2, rc);
    lambda3(idx) = cast(l3, rc);
end
clear Hxx Hxy Hxz Hyy Hyz Hzz idx l2 l3

lambda2(~isfinite(lambda2) | abs(lambda2) < eigTol) = 0;
lambda3(~isfinite(lambda3) | abs(lambda3) < eigTol) = 0;

% ---- (3) 3-factor regularisation of lambda3 ----
tube = lambda2 < 0 & lambda3 < 0;
m3 = cast(0, rc);
if any(tube(:))
    m3 = min(lambda3(tube));
end

lambda3t = lambda3;
lambda4  = lambda3;
if m3 < 0
    lambda3t(lambda3 < 0 & lambda3 >= opts.tau*m3)  = opts.tau*m3;
    lambda4 (lambda3 < 0 & lambda3 >= opts.tau2*m3) = opts.tau2*m3;
end

% ---- (4) FA-like measure ----
a2 = abs(lambda2);  a3 = abs(lambda3t);  a4 = abs(lambda4);
mu = (a2 + a3 + a4) ./ cast(3, rc);
numer = sqrt((a2 - mu).^2 + (a3 - mu).^2 + (a4 - mu).^2);
denom = sqrt(a2.^2 + a3.^2 + a4.^2);
denom(denom == 0) = eps0;
fa = numer ./ denom;
fa(~isfinite(fa)) = 0;

% ---- pack ----
geom.lambda2  = lambda2;
geom.lambda3  = lambda3;
geom.lambda3t = lambda3t;
geom.lambda4  = lambda4;
geom.fa       = fa;
end
