function R = vesselnessResponse(Lambda1, Lambda2, beta, c)
%VESSELNESSRESPONSE Frangi vesselness (2D)
%
%   R = vesselnessResponse(Lambda1, Lambda2, BETA, C)
%
% DESCRIPTION
%   Computes the classic Frangi vesselness measure for 2D images
%   using Hessian eigenvalues.
%
%   This implementation enforces the Frangi eigenvalue convention
%   locally (|Lambda1| <= |Lambda2|) and is robust to upstream eigen ordering.
%
% PARAMETERS
%   BETA  - controls blob vs tube discrimination (typical: 0.3–0.7)
%   C     - controls noise suppression (typical: 10–20)
%
% GUARANTEES
%   - strong response on tubular structures
%   - suppression of blobs relative to tubes
%   - suppression of flat background
%
% NON-GOALS
%   - thin > thick discrimination
%   - binary segmentation
%
% REFERENCE
%   Frangi et al., "Multiscale Vessel Enhancement Filtering",
%   MICCAI, 1998.


% Ensure |Lambda1| <= |Lambda2|
swap = abs(Lambda1) > abs(Lambda2);
Lambda1s = Lambda1;
Lambda2s = Lambda2;
Lambda1s(swap) = Lambda2(swap);
Lambda2s(swap) = Lambda1(swap);

% Avoid division by zero
Lambda1s(Lambda1s == 0) = eps;

% Frangi measures
Rb = (Lambda1s ./ Lambda2s).^2;      % NOTE: small over large
S2 = Lambda1s.^2 + Lambda2s.^2;

R = exp(-Rb / (2*beta^2)) .* (1 - exp(-S2 / (2*c^2)));

% Polarity: bright vessels on dark background
R(Lambda2s >= 0) = 0;
end
