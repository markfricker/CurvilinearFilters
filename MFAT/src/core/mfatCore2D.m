function geom = mfatCore2D(I, sigma, opts)
% =========================================================================
% MFAT – Multiscale Fractional Anisotropy Tensor Framework
%
% Author:
%   H. Alhasson, M. Alharbi, B. Obara
%
% Refactored framework & extensions:
%   MD Fricker, Jan 2026
%
% Citation:
%   H. Alhasson, M. Alharbi, B. Obara,
%   "2D and 3D Vascular Structures Enhancement via
%    Multiscale Fractional Anisotropy Tensor",
%   ECCV Workshops (BioImage Computing), 2018.
%
% Background:
%   P. J. Basser et al.,
%   "MR Diffusion Tensor Spectroscopy and Imaging",
%   Biophysical Journal, 1994.
%
% Related work:
%   A. F. Frangi et al.,
%   "Multiscale Vessel Enhancement Filtering",
%   MICCAI, 1998.
%
% License:
%   Academic / research use. Please cite the above work.
% =========================================================================

% MFATCORE2D  MFAT geometry core (single scale)
%
% OVERVIEW
%   Computes MFAT geometric quantities at a single scale:
%     - Hessian principal eigenvalue,
%     - three-factor eigenvalue regularisation,
%     - fractional-anisotropy-like tensor measure.
%
%   This function performs NO thresholding, NO aggregation, and NO
%   decision-making. It is purely geometric.
%
% INPUTS
%   I     - 2D image (already normalized).
%   sigma - Gaussian scale.
%   opts  - Struct with fields:
%           .tau, .tau2, .whiteOnDark, .precision
%
% OUTPUT
%   geom - Struct containing:
%           .lambda2, .lambda3, .lambda4, .fa
%
% DESIGN NOTE
%   This separation ensures that MFAT-λ, MFAT-Prob, entropy, and fractional
%   variants all operate on the same geometric backbone.
%
% REFERENCES
%   [1] H. Alhasson et al., ECCV Workshops, 2018.


% ---- defaults ----
if ~isfield(opts,'whiteOnDark'), opts.whiteOnDark = true; end
if ~isfield(opts,'tau'),         opts.tau = 0.03;        end
if ~isfield(opts,'tau2'),        opts.tau2 = 0.3;        end
if ~isfield(opts,'precision'),   opts.precision = 'single'; end

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

I = cast(I, rc);

% ---- (1) Hessian eigenvalue ----
[~, lambda2] = imageEigenvalues2D( ...
    I, sigma, opts.whiteOnDark, rc, eigTol);

% ---- (2) 3-factor regularisation ----
lambda3 = lambda2;
m3 = min(lambda3(:));
if m3 < 0
    lambda3(lambda3 < 0 & lambda3 >= opts.tau*m3) = opts.tau*m3;
end

lambda4 = lambda2;
m4 = min(lambda4(:));
if m4 < 0
    lambda4(lambda4 < 0 & lambda4 >= opts.tau2*m4) = opts.tau2*m4;
end

% ---- (3) FA-like response ----
mu = (abs(lambda2) + abs(lambda3) + abs(lambda4)) ./ cast(3, class(lambda2));


numer = sqrt( ...
    (abs(lambda2) - mu).^2 + ...
    (abs(lambda3) - mu).^2 + ...
    (abs(lambda4) - mu).^2 );

denom = sqrt(lambda2.^2 + lambda3.^2 + lambda4.^2);
denom(denom == 0) = eps0;

fa = numer ./ denom;
fa(~isfinite(fa)) = 0;

% ---- pack ----
geom.lambda2 = lambda2;
geom.lambda3 = lambda3;
geom.lambda4 = lambda4;
geom.fa      = fa;
end
