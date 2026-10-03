function vesselness = mfatResponseLambda3D(geom, vesselness, opts)
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

% MFATRESPONSELAMBDA3D  MFAT-λ response update for a volume (single scale)
%
% OVERVIEW
%   3D counterpart of mfatResponseLambda: response = 1 - sqrt(3/2)*FA,
%   kept only where lambda2 < 0 and lambda3 < 0 (a bright tube: two
%   strongly negative cross-sectional curvatures), then the same D·tanh +
%   max multiscale fusion as 2D.
%
%   The published 3D rule block (FractionalIstropicTensor3D.m) reduces to
%   exactly this: with x = lambda3t - lambda2, "lambda3t > x" is
%   equivalent to lambda2 > 0, and "x == max(x)" only ever selects voxels
%   the lambda2 >= 0 rule zeroes again. Those two lines are therefore
%   omitted. Note that the 3D rules differ from the 2D ones (2D also keeps
%   only voxels past the tau clip): this follows the published code.
%
%   A sheet (|lambda2| ~ 0, |lambda3| large) gets sqrt(3/2)*FA ~ 0.71,
%   i.e. a response of ~0.29 at most; a round tube with |lambda2| ~ |lambda3| beyond the
%   tau2 clip gets ~1. tau/tau2 thus set a soft contrast ramp relative to
%   the strongest tube in the volume.
%
% INPUTS
%   geom       - Geometry struct from mfatCore3D.
%   vesselness - Accumulated response (empty on first scale).
%   opts       - Struct with fields .D and .precision.
%
% OUTPUT
%   vesselness - Updated response in [0,1].
%
% See also: mfatResponseLambda, mfatCore3D

scale = sqrt(cast(3, class(geom.fa)) / cast(2, class(geom.fa)));
response = imcomplementSafe(scale .* geom.fa, opts.precision);

response(geom.lambda2 >= 0 | geom.lambda3 >= 0) = 0;
response(~isfinite(response)) = 0;

if isempty(vesselness)
    vesselness = response;
else
    vesselness = vesselness + opts.D .* tanh(response - opts.D);
    vesselness = max(vesselness, response);
end

vesselness = min(max(vesselness,0),1);
end
