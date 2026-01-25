function R = ridgeResponse(l1,l2,alpha)
%RIDGERESPONSE Ridge detector (Sato-style)
%
%   R = ridgeResponse(L1, L2, ALPHA)
%
% DESCRIPTION
%   Computes a ridge response based on curvature magnitude and
%   anisotropy of the Hessian.
%
% PARAMETERS
%   ALPHA - anisotropy penalty (typical: ~0.5)
%
% GUARANTEES
%   - responds to elongated structures
%   - suppresses flat background
%
% NON-GOALS
%   - blob suppression
%   - thin > thick preference
%
% REFERENCE
%   Sato et al., "Three-dimensional multi-scale line filter",
%   Medical Image Analysis, 1998.


R = abs(l2) .* exp(-(l1./(abs(l2)+eps)).^2/(2*alpha^2));
end
