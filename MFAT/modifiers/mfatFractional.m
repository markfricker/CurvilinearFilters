function R_frac = mfatFractional(R_lambda, varargin)
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

% MFATFRACTIONAL  Fractional shaping of MFAT-λ response
%
% OVERVIEW
%   Applies a fractional power-law nonlinearity to MFAT-λ to reshape
%   response contrast without altering geometry.
%
% PARAMETER
%   'alpha' - Fractional exponent.
%             alpha < 1 boosts weak responses.
%             alpha > 1 sharpens strong responses.


p = inputParser;
addParameter(p,'alpha',0.7,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:});

alpha = p.Results.alpha;

R_frac = R_lambda .^ alpha;

mx = max(R_frac(:));
if mx > 0
    R_frac = R_frac ./ mx;
end

R_frac = min(max(R_frac,0),1);
end
