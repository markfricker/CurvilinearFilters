function out = mfatProb(I, sigmas, varargin)
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

% MFATPROB  Probabilistic MFAT based on MFAT-λ
%
% OVERVIEW
%   MFAT-Prob wraps the deterministic MFAT-λ response in a probabilistic
%   framework by modelling per-scale responses with Beta distributions
%   and accumulating log-likelihood ratios.
%
%   It does NOT redefine the detector; it provides calibrated confidence.
%
% OUTPUT
%   out - Posterior probability of curvilinear structure.
%
% REFERENCES
%   [1] H. Alhasson et al., ECCV Workshops, 2018.

% ---- parse ----
p = inputParser;
addParameter(p,'tau',0.03);
addParameter(p,'tau2',0.3);
addParameter(p,'D',0.27);
addParameter(p,'whiteOnDark',true);
addParameter(p,'precision','single');
parse(p,varargin{:});
opts = p.Results;

% ---- precision ----
if strcmpi(opts.precision,'single')
    eps0 = eps('single');
else
    eps0 = eps('double');
end

% ---- preprocess ----
I = cast(I, opts.precision);
I = I ./ (max(I(:)) + eps0);
sigmas = double(sigmas(:)');

% ---- initialise probabilistic state ----
state = mfatResponseProb('init', size(I), [], opts);

% ---- multiscale loop ----
for k = 1:numel(sigmas)
    % MFAT core
    geom = mfatCore2D(I, sigmas(k), opts);

    % MFAT-λ response at this scale
    lambdaResp = mfatResponseLambda(geom, [], opts);

    % accumulate probabilistic evidence
    state = mfatResponseProb('accumulate', lambdaResp, state, opts);
end

% ---- final posterior ----
out = mfatResponseProb('finalize', [], state, opts);

% optional normalization (already in [0,1])
out = min(max(out,0),1);
end
