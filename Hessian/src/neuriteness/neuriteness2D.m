function [response, direction, Vx, Vy] = neuriteness2D(I, sigma, Precision)
%NEURITENESS2D  Meijering neuriteness filter (2D, single scale)
%
%   [R, DIR] = neuriteness2D(I, SIGMA)
%   [R, DIR, VX, VY] = neuriteness2D(...)
%
% DESCRIPTION
%   Implements the 2D neuriteness filter described in:
%
%     E. Meijering et al.,
%     "Design and Validation of a Tool for Neurite Tracing and Analysis
%      in Fluorescence Microscopy Images",
%     Cytometry Part A, 58A:167–176, 2004
%
%   This is a SINGLE-SCALE filter. No multiscale aggregation is performed.
%
%   The algorithm:
%     1) computes the Hessian matrix at scale SIGMA;
%     2) performs Lindeberg scale normalization;
%     3) applies the Meijering modified-Hessian transformation;
%     4) computes neuriteness from the dominant eigenvalue;
%     5) normalizes responses globally.
%
% INPUTS
%   I         - 2D grayscale image
%   SIGMA     - Gaussian scale (scalar)
%   PRECISION - 'single' or 'double' (default: 'single')
%
% OUTPUTS
%   R         - neuriteness response image (range [0,1])
%   DIR       - local neurite orientation (radians, tangent direction)
%   VX, VY    - x/y components of the tangent direction vector
%
% NOTES
%   - Neuriteness is a pure shape measure (no contrast weighting).
%   - Bright neurites on dark background are assumed.
%   - Polarity suppression is intrinsic (L1 < 0).
%   - For multiscale detection, combine with vesselness explicitly.
%
% SEE ALSO
%   hessian2DFilters, Hessian2D, eig2image

if nargin < 3
    Precision = 'single';
end

% -------------------------------------------------------------------------
% Precision handling
% -------------------------------------------------------------------------
if strcmpi(Precision,'single')
    I = single(I);
else
    I = double(I);
end

% -------------------------------------------------------------------------
% Hessian eigenvalues and eigenvectors (shared core)
% -------------------------------------------------------------------------
[L1, L2, Vnx, Vny] = hessianEigen2D(I, sigma, Precision);

% -------------------------------------------------------------------------
% Modified Hessian (Meijering)
% -------------------------------------------------------------------------
alpha = -1/3;

L1t = L1;
L2t = L2;

L1 = L1t + alpha * L2t;
L2 = L2t + alpha * L1t;

% -------------------------------------------------------------------------
% Sort so |L1| >= |L2| (dominant eigenvalue first)
% -------------------------------------------------------------------------
idx = abs(L1) < abs(L2);

L1s = L1;
L2s = L2;

L1s(idx) = L2(idx);
L2s(idx) = L1(idx);

% -------------------------------------------------------------------------
% Neuriteness response (global normalization)
% -------------------------------------------------------------------------
Lmin = min(L1s(:));

response = zeros(size(I), 'like', I);
mask = (L1s < 0);

response(mask) = L1s(mask) / Lmin;

% -------------------------------------------------------------------------
% Direction: neurite tangent
% -------------------------------------------------------------------------
% Hessian eigenvector (Vnx, Vny) is normal to the structure.
% Tangent direction is perpendicular to it.

Vx = -Vny;
Vy =  Vnx;

direction = atan2(Vy, Vx);

end
