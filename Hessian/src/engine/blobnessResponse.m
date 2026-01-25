function R = blobnessResponse(l1,l2)
%BLOBNESSRESPONSE Blob detector (Determinant-of-Hessian-like)
%
%   R = blobnessResponse(L1, L2)
%
% DESCRIPTION
%   Computes a blobness response based on the determinant of the
%   Hessian matrix. This response enhances isotropic structures
%   where curvature is similar in all directions.
%
%   The response is high at blob centers and low in flat regions.
%
% INPUTS
%   L1, L2 - Hessian eigenvalues at a given scale
%
% OUTPUT
%   R      - blobness response (non-negative)
%
% GUARANTEES
%   - strong response at blob centers
%   - suppresses flat background regions
%
% DOES NOT GUARANTEE
%   - zero response on lines or ridges
%   - discrimination between blobs and thick lines
%
% NOTES
%   Blobness is best interpreted comparatively across scales.
%   Use a sigma range that spans expected blob radii.
%
% REFERENCE
%   T. Lindeberg, "Feature detection with automatic scale selection",
%   International Journal of Computer Vision, 1998.


R = abs(l1 .* l2);
end
