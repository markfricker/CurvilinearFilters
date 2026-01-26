function vesselness = mfatResponseLambda(geom, vesselness, opts)
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

% MFATRESPONSELAMBDA  Deterministic MFAT-λ response update (single scale)
%
% OVERVIEW
%   Applies strict MFAT post-processing rules and deterministic multiscale
%   aggregation to geometric quantities computed by mfatCore2D.
%
% INPUTS
%   geom       - Geometry struct from mfatCore2D.
%   vesselness - Accumulated MFAT response (empty on first scale).
%   opts       - Struct with field .D and .precision.
%
% OUTPUT
%   vesselness - Updated MFAT-λ response.
%
% NOTES
%   The logic in this function reproduces the exact behaviour described in
%   Alhasson et al. and should not be altered lightly.


% ---- precision ----
if strcmpi(opts.precision,'single')
    eps0 = eps('single');
else
    eps0 = eps('double');
end

lambda2 = geom.lambda2;
lambda3 = geom.lambda3;
fa      = geom.fa;

% ---- inversion / scaling ----
scale = sqrt(cast(3, class(fa)) / cast(2, class(fa)));
response = imcomplementSafe(scale .* fa, opts.precision);


% ---- strict MFAT rules ----
x = lambda3 - lambda2;
response(x == min(x(:))) = 1;
response(x < max(x(:)))  = 0;
response(lambda2 >= 0)   = 0;
response(lambda3 >= 0)   = 0;
response(~isfinite(response)) = 0;

% ---- multiscale aggregation ----
if isempty(vesselness)
    vesselness = response;
else
    vesselness = vesselness + opts.D .* tanh(response - opts.D);
    vesselness = max(vesselness, response);
end

vesselness = min(max(vesselness,0),1);
end
