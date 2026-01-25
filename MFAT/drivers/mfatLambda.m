function out = mfatLambda(I, sigmas, varargin)
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

% MFATLAMBDA  Deterministic Multiscale Fractional Anisotropy Tensor (MFAT-λ)
%
%   out = mfatLambda(I, sigmas)
%   out = mfatLambda(I, sigmas, 'Name', Value, ...)
%
% OVERVIEW
%   MFAT-λ implements the deterministic Multiscale Fractional Anisotropy
%   Tensor filter proposed by Alhasson et al. for enhancing curvilinear
%   structures. The method combines:
%     - Hessian-based curvature analysis,
%     - three-factor eigenvalue regularisation,
%     - fractional-anisotropy-like tensor measures,
%     - strict post-processing rules,
%     - deterministic multiscale aggregation.
%
%   The output is a near-binary curvilinear evidence map.
%
% ALGORITHM
%   For each scale σ in sigmas:
%     1. Compute Hessian eigenvalues at scale σ.
%     2. Form a three-factor tensor using tau and tau2 regularisation.
%     3. Compute a fractional-anisotropy-like response.
%     4. Apply strict MFAT post-processing rules.
%     5. Aggregate responses across scales using D·tanh and max fusion.
%
% INPUTS
%   I       - 2D grayscale image (numeric).
%   sigmas  - Vector of Gaussian scales.
%
% NAME-VALUE PARAMETERS
%   'tau'         - Lower eigenvalue clipping factor (default 0.03).
%   'tau2'        - Upper eigenvalue clipping factor (default 0.3).
%   'D'           - Multiscale aggregation strength (default 0.27).
%   'whiteOnDark' - true if structures are bright on dark background.
%   'precision'   - 'single' (default) or 'double'.
%
% OUTPUT
%   out - MFAT-λ response in [0,1].
%
% PARAMETER TUNING
%   - Thin structures: smaller tau (0.005–0.03), smaller sigmas.
%   - Thicker structures: larger sigmas, larger tau2.
%   - Crisper output: increase D.
%
% REFERENCES
%   [1] H. Alhasson et al., ECCV Workshops, 2018.
%   [2] P. J. Basser et al., Biophysical Journal, 1994.
%   [3] A. F. Frangi et al., MICCAI, 1998.

% ---- parse ----
p = inputParser;
addParameter(p,'tau',0.03);
addParameter(p,'tau2',0.3);
addParameter(p,'D',0.27);
addParameter(p,'whiteOnDark',true);
addParameter(p,'precision','single');
parse(p,varargin{:});
opts = p.Results;

% ---- precision-safe eps ----
if strcmpi(opts.precision,'single')
    eps0 = eps('single');
else
    eps0 = eps('double');
end

% ---- preprocess ----
I = cast(I, opts.precision);
I = I ./ (max(I(:)) + eps0);
sigmas = double(sigmas(:)');

vesselness = [];

% ---- multiscale MFAT ----
for k = 1:numel(sigmas)
    geom = mfatCore2D(I, sigmas(k), opts);
    vesselness = mfatResponseLambda(geom, vesselness, opts);
end

% ---- final cleanup ----
out = vesselness;
out = out ./ (max(out(:)) + eps0);
out(out < cast(1e-2, class(out))) = cast(0, class(out));
end
