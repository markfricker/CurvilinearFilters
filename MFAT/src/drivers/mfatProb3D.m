function out = mfatProb3D(V, sigmas, varargin)
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

% MFATPROB3D  Probabilistic MFAT for a 3D volume
%
%   out = mfatProb3D(V, sigmas, 'Name', Value, ...)
%
% OVERVIEW
%   3D counterpart of mfatProb: each scale's MFAT-λ response
%   (mfatCore3D + mfatResponseLambda3D) is scored under Beta
%   vessel/background models and the log-likelihood ratios are summed
%   across scales (mfatResponseProb, shared unchanged with 2D), then
%   mapped to a posterior.
%
%   As in 2D, this is the framework's Beta-LLR wrapper, NOT the published
%   "PFAT" formula (ProbabiliticFractionalIstropicTensor3D.m), which
%   instead normalises the eigenvalues by their trace before the FA step.
%   'D' is accepted for parameter compatibility with mfatLambda3D but not
%   used (scales are fused by LLR summation, not D·tanh).
%
% INPUTS / NAME-VALUE PARAMETERS
%   As mfatLambda3D.
%
% OUTPUT
%   out - Posterior probability of tubular structure, in [0,1].
%
% See also: mfatProb, mfatLambda3D, mfatResponseProb

% ---- parse ----
p = inputParser;
addParameter(p,'tau',0.03);
addParameter(p,'tau2',0.3);
addParameter(p,'D',0.27);
addParameter(p,'whiteOnDark',true);
addParameter(p,'precision','single');
addParameter(p,'spacing',[1 1 1]);
parse(p,varargin{:});
opts = p.Results;

if strcmpi(opts.precision,'single')
    eps0 = eps('single');
else
    eps0 = eps('double');
end

% ---- preprocess ----
V = cast(V, opts.precision);
V = V ./ (max(V(:)) + eps0);
if iscell(sigmas)
    sigmas = sigmas(:)';
else
    sigmas = num2cell(double(sigmas(:)'));
end

state = mfatResponseProb('init', size(V), [], opts);

for k = 1:numel(sigmas)
    geom = mfatCore3D(V, sigmas{k}, opts);
    lambdaResp = mfatResponseLambda3D(geom, [], opts);
    state = mfatResponseProb('accumulate', lambdaResp, state, opts);
end

out = mfatResponseProb('finalize', [], state, opts);
out = min(max(out,0),1);
end
