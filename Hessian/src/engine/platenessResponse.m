function R = platenessResponse(l1,l2,alpha)
%PLATENESSRESPONSE Plateness filter (2D analogue)
%
%   R = platenessResponse(L1, L2, ALPHA)
%
% DESCRIPTION
%   Computes a plateness response for 2D images, designed to enhance
%   thick elongated regions that are intermediate between ridges and
%   blobs.
%
%   In 2D, plateness is a heuristic analogue of sheetness in 3D and
%   is most useful for membranes, thick bands, or wide elongated
%   structures.
%
% INPUTS
%   L1, L2 - Hessian eigenvalues at a given scale
%   ALPHA  - anisotropy tolerance parameter (typical: ~0.5)
%
% OUTPUT
%   R      - plateness response (non-negative)
%
% GUARANTEES
%   - higher response on thick elongated structures than thin ones
%   - suppresses flat background regions
%
% DOES NOT GUARANTEE
%   - blob suppression
%   - thin-line enhancement
%
% NOTES
%   Plateness is weakly defined in 2D and should be interpreted
%   qualitatively rather than as a strict classifier.
%
% REFERENCE
%   Y. Sato et al., "Three-dimensional multi-scale line filter",
%   Medical Image Analysis, 1998.


R = abs(l2) .* (1 - exp(-(abs(l1)./(abs(l2)+eps)).^2/(2*alpha^2)));
end
