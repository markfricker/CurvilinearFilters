function R = platenessResponse3D(L1, L2, L3, alpha, beta, c)
%PLATENESSRESPONSE3D  Sheetness (3D, planar/membrane structures)
%
%   R = platenessResponse3D(L1, L2, L3, ALPHA, BETA, C)
%
% DESCRIPTION
%   Genuine 3D sheet/plate detection -- unlike the 2D engine's
%   platenessResponse.m, which is explicitly documented there as "weakly
%   defined in 2D, interpret qualitatively" (a heuristic stand-in for
%   "thick region", not real sheetness), this measure has an actual
%   geometric basis: "flatness" requires a third dimension to be flat
%   *in*, so a real sheet has TWO small eigenvalues (the in-plane
%   directions, low curvature along the sheet) and ONE large eigenvalue
%   (across the sheet thickness). A 2D Hessian, with only two eigenvalues,
%   cannot express this at all -- this is the reason a native 3D Enhance
%   stage is needed for ER cisternae, not an optimisation on top of an
%   already-adequate 2D pipeline.
%
%   Uses the SAME Ra/Rb/S primitives as vesselnessResponse3D -- a tube has
%   Ra=|L2|/|L3| near 1 (L2 and L3 comparable: two directions of curvature
%   across the tube), a sheet has Ra near 0 (L2 much smaller than L3: only
%   one dominant direction of curvature, across the sheet). Plateness is
%   the same construction as vesselness with that one factor inverted --
%   the way Frangi-style measures have been adapted for sheetness in the
%   literature (Descoteaux et al. 2006 built a multi-scale Hessian sheet
%   measure for thin bone in CT along these lines; their blob term differs
%   in detail from the Rb used here).
%
%   Requires |L1| <= |L2| <= |L3| (eig3volume's convention).
%
%   R = exp(-Ra^2/2*alpha^2) * exp(-Rb^2/2*beta^2) * (1 - exp(-S^2/2*c^2))
%
%   Polarity is NOT applied here -- see vesselnessResponse3D's header;
%   the caller (hessian3DFilters) applies the WhiteOnDark sign gate.
%
% PARAMETERS
%   ALPHA, BETA, C -- same roles and typical values as vesselnessResponse3D
%
% REFERENCE
%   Frangi A.F. et al. (1998) MICCAI 1998, LNCS 1496:130-137 (Ra/Rb/S
%   construction). Descoteaux M., Audette M., Chinzei K., Siddiqi K. (2006)
%   "Bone enhancement filtering: application to sinus bone segmentation and
%   simulation of pituitary surgery", Computer Aided Surgery 11(5):247-255,
%   doi:10.3109/10929080601017212 (Hessian sheetness measure; conference
%   version MICCAI 2005, LNCS 3749:9-16). [Corrected 2026-10-01: this
%   reference previously gave a vessel-segmentation title with Medical Image
%   Analysis 10(6):850-862, which is an unrelated paper.]
%
% See also: vesselnessResponse3D, hessian3DFilters, eig3volume

L2safe = L2; L2safe(L2safe == 0) = eps;
L3safe = L3; L3safe(L3safe == 0) = eps;

Ra = abs(L2) ./ abs(L3safe);
Rb = abs(L1) ./ sqrt(abs(L2safe .* L3safe));
S  = sqrt(L1.^2 + L2.^2 + L3.^2);

R = exp(-(Ra.^2) / (2*alpha^2)) .* ...
    exp(-(Rb.^2) / (2*beta^2)) .* ...
    (1 - exp(-(S.^2) / (2*c^2)));
end
