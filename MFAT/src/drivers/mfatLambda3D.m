function out = mfatLambda3D(V, sigmas, varargin)
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

% MFATLAMBDA3D  Deterministic MFAT-λ for a 3D volume
%
%   out = mfatLambda3D(V, sigmas)
%   out = mfatLambda3D(V, sigmas, 'Name', Value, ...)
%
% OVERVIEW
%   3D counterpart of mfatLambda, on the native (possibly anisotropic)
%   voxel grid. For each scale: Hessian eigenvalues (mfatCore3D), the
%   three-factor FA response (mfatResponseLambda3D), D·tanh + max fusion.
%   See mfatCore3D for how this differs from the published 3D code.
%
% INPUTS
%   V       - 3D volume (numeric), [dim1 dim2 dim3].
%   sigmas  - Scales in PHYSICAL units (same units as 'spacing'):
%             numeric vector (one isotropic physical scale per step), or a
%             cell array of [s1 s2 s3] vectors (independent per-axis scale
%             per step -- needed when Z resolution is much coarser than
%             XY; see applyHessian3DAniso). With the default spacing
%             [1 1 1] these are simply pixel sigmas.
%
% NAME-VALUE PARAMETERS
%   'tau'         - Lower eigenvalue clipping factor (default 0.03).
%   'tau2'        - Upper eigenvalue clipping factor (default 0.3).
%   'D'           - Multiscale aggregation strength (default 0.27).
%   'whiteOnDark' - true if structures are bright on dark (default true).
%   'precision'   - 'single' (default) or 'double'.
%   'spacing'     - Voxel size [d1 d2 d3] (default [1 1 1]).
%
% OUTPUT
%   out - MFAT-λ response in [0,1], same size as V.
%
% See also: mfatLambda, mfatProb3D, mfatCore3D

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
sigmas = localSigmaCell(sigmas);

vesselness = [];

% ---- multiscale MFAT ----
for k = 1:numel(sigmas)
    geom = mfatCore3D(V, sigmas{k}, opts);
    vesselness = mfatResponseLambda3D(geom, vesselness, opts);
end

% ---- final cleanup ----
out = vesselness;
out = out ./ (max(out(:)) + eps0);
out(out < cast(1e-2, class(out))) = cast(0, class(out));
end

function c = localSigmaCell(sigmas)
% numeric vector -> one cell per scale step; cell array passed through
if iscell(sigmas)
    c = sigmas(:)';
else
    c = num2cell(double(sigmas(:)'));
end
end
