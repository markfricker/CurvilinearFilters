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
%   Descoteaux et al. (2006) sheetness. Ra=|L2|/|L3| is the same primitive
%   as vesselnessResponse3D -- near 1 for a tube (two directions of
%   curvature across it), near 0 for a sheet (one, across its thickness)
%   -- with the factor inverted. The blob term is Descoteaux's
%     Rb = |2|L3| - |L2| - |L1|| / |L3|     (sheet ~2, tube ~1, blob ~0)
%   used as (1 - exp(...)), so a blob is driven to zero.
%   [Changed 2026-10-03 from the vesselness-style |L1|/sqrt(|L2*L3|),
%   exp(-Rb^2...): that ratio is bounded (Rb^2 <= Ra, since |L1|<=|L2|),
%   so it never destabilised, but on a sheet it is a ratio of two noise
%   eigenvalues and barely suppressed blobs (term 0.20 vs 0.02 now, on a
%   blurred-sphere phantom; sheet term 0.96 -> 1.00, tube 1.00 -> 0.89).]
%
%   Requires |L1| <= |L2| <= |L3| (eig3volume's convention).
%
%   R = exp(-Ra^2/2*alpha^2) * (1 - exp(-Rb^2/2*beta^2)) * (1 - exp(-S^2/2*c^2))
%
%   Polarity is NOT applied here -- see vesselnessResponse3D's header;
%   the caller (hessian3DFilters) applies the WhiteOnDark sign gate, on L3
%   only for this measure (the in-plane L1/L2 have a noise-driven sign).
%
% PARAMETERS
%   ALPHA, BETA, C -- same roles and typical values as vesselnessResponse3D
%   (Descoteaux use alpha = beta = 0.5; c data-dependent)
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

L3safe = L3; L3safe(L3safe == 0) = eps;

Ra = abs(L2) ./ abs(L3safe);
Rb = abs(2*abs(L3) - abs(L2) - abs(L1)) ./ abs(L3safe);
S  = sqrt(L1.^2 + L2.^2 + L3.^2);

R = exp(-(Ra.^2) / (2*alpha^2)) .* ...
    (1 - exp(-(Rb.^2) / (2*beta^2))) .* ...
    (1 - exp(-(S.^2) / (2*c^2)));
end
