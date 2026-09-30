function R = vesselnessResponse3D(L1, L2, L3, alpha, beta, c)
%VESSELNESSRESPONSE3D  Frangi vesselness (3D, tubular structures)
%
%   R = vesselnessResponse3D(L1, L2, L3, ALPHA, BETA, C)
%
% DESCRIPTION
%   The genuine three-eigenvalue Frangi 1998 vesselness measure. Unlike
%   the 2D engine's vesselnessResponse.m (which only has two eigenvalues
%   and so cannot distinguish tube from plate at all), this is the actual
%   formula Frangi's paper was about -- 2D "vesselness" is really just
%   "not blob", whereas 3D vesselness genuinely discriminates tube from
%   plate from blob using all three eigenvalues.
%
%   Requires |L1| <= |L2| <= |L3| (eig3volume's convention).
%
%   Ra = |L2|/|L3|              -- tube (~1) vs plate (~0) discriminator
%   Rb = |L1|/sqrt(|L2*L3|)     -- blob (large) vs line/plate (~0)
%   S  = sqrt(L1^2+L2^2+L3^2)   -- structureness (background suppression)
%
%   R = (1 - exp(-Ra^2/2*alpha^2)) * exp(-Rb^2/2*beta^2) * (1 - exp(-S^2/2*c^2))
%
%   Polarity (WhiteOnDark: does L2/L3<0 or >0 mean "bright structure") is
%   NOT applied here -- R is polarity-symmetric (depends on L2/L3 only
%   through ratios/squares), exactly like the 2D engine's response
%   functions. The caller (hessian3DFilters) applies the sign gate, so one
%   pass of eigenvalues serves either polarity without recomputing them.
%
% PARAMETERS
%   ALPHA - tube vs plate discrimination (typical: ~0.5, Frangi's default)
%   BETA  - blob suppression (typical: ~0.5)
%   C     - noise/background suppression (typical: half the max Frobenius
%           norm seen in the image; same role as the 2D engine's C)
%
% REFERENCE
%   Frangi A.F. et al. (1998) "Multiscale Vessel Enhancement Filtering",
%   MICCAI 1998, LNCS 1496:130-137. https://doi.org/10.1007/BFb0056195
%
% See also: platenessResponse3D, hessian3DFilters, eig3volume

L2safe = L2; L2safe(L2safe == 0) = eps;
L3safe = L3; L3safe(L3safe == 0) = eps;

Ra = abs(L2) ./ abs(L3safe);
Rb = abs(L1) ./ sqrt(abs(L2safe .* L3safe));
S  = sqrt(L1.^2 + L2.^2 + L3.^2);

R = (1 - exp(-(Ra.^2) / (2*alpha^2))) .* ...
    exp(-(Rb.^2) / (2*beta^2)) .* ...
    (1 - exp(-(S.^2) / (2*c^2)));
end
