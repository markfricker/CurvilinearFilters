function R_ent = mfatEntropyWeight(R_lambda, faStack, varargin)
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

% MFATENTROPYWEIGHT  Entropy-based weighting of MFAT-λ
%
% OVERVIEW
%   Applies a conservative entropy-based weighting to MFAT-λ responses,
%   downweighting regions of ambiguous tensor anisotropy.
%
% INPUTS
%   R_lambda - MFAT-λ response.
%   faStack  - FA maps per scale (H×W×S).
%
% PARAMETERS
%   'beta'    - Entropy strength (default 0.5).
%   'combine' - Aggregation across scales: 'mean', 'min', 'max'.
%
% USE CASES
%   - Suppress spurious junctions.
%   - Improve reliability ranking.


p = inputParser;
addParameter(p,'beta',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'combine','mean',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
beta = p.Results.beta;
combine = p.Results.combine;

% ---- stack handling ----
if iscell(faStack)
    faAll = cat(3, faStack{:});
else
    faAll = faStack;
end

% ---- entropy per scale ----
eps0 = eps(class(faAll));
H = -faAll .* log(faAll + eps0);   % entropy surrogate

% ---- aggregate entropy across scales ----
switch lower(combine)
    case 'mean'
        Hagg = mean(H,3);
    case 'min'
        Hagg = min(H,[],3);
    case 'max'
        Hagg = max(H,[],3);
    otherwise
        error('Unknown combine mode.');
end

% ---- entropy weight ----
wEnt = exp(-beta * Hagg);

% ---- apply conservatively ----
R_ent = R_lambda .* wEnt;

% clamp
R_ent = min(max(R_ent,0),1);
end
