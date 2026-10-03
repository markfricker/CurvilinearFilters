function [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3DAniso(I, sigmaPhysical, spacing)
%APPLYHESSIAN3DANISO  3D Hessian on an anisotropic voxel grid, in physical units.
%
%   [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3DAniso(I, sigmaPhysical, spacing)
%
% Computes the Hessian directly on the NATIVE (anisotropic) voxel grid --
% no resampling -- by (1) smoothing with a per-axis PIXEL sigma chosen so
% the smoothing scale matches sigmaPhysical in every direction despite the
% unequal voxel spacing, then (2) rescaling the resulting derivatives into
% PHYSICAL curvature units (per micron^2, or whatever unit `spacing` is
% in) so entries computed along different array axes are directly
% comparable -- the spacing-aware correction ITK's own Hessian filters
% apply for anisotropic input, not a novel trick.
%
% This is the alternative to isotropic-resampling the greyscale volume
% before calling applyHessian3D/hessianEigen3D: on a real 76-slice, 5.3x
% anisotropic ER volume, isotropic resampling inflated Z to 404 slices and
% cost ~450s per multiscale vesselness pass; this function works on the
% native 76-slice grid directly.
%
% PER-AXIS SIGMA (added after real-data testing): sigmaPhysical may be a
% SCALAR (one physical smoothing scale applied through all three axes --
% the original behaviour) OR a 3-element vector [sigma1 sigma2 sigma3]
% giving an INDEPENDENT physical scale per axis. This matters when lateral
% and axial resolution are too different for one scale to serve both: on
% real 5.3x-anisotropic ER data, forcing enough Z-sigma to get even 1
% non-degenerate Z-pixel of smoothing (sigmaPhysical >= dz) forced the SAME
% scalar through XY too, requiring >5 pixels of XY smoothing there --
% coarser than the tubules themselves, so the filter stopped resolving
% fine tubules and started responding to much larger rounded structures
% instead (visually confirmed: vesselness went from tracing thin tubules
% to tracing blob/crescent rims once s3 was forced >=1 with one shared
% sigma). A per-axis vector lets XY stay matched to the true tubule width
% while Z independently uses whatever is actually resolvable.
%
% AXIS CONVENTION -- matches applyHessian3D.m exactly: "Dxx" is the 2nd
% derivative along array dimension 1, "Dyy" along dimension 2, "Dzz" along
% dimension 3 (an internal x/y/z labelling choice that dates from the 2D
% code's row/column ndgrid usage -- it does not need to mean anything
% about real-world X/Y/Z, since eig3volume and the response functions are
% rotation-invariant and never interpret which axis is which). `spacing`
% follows the SAME dimension order: spacing(1) is dim 1's voxel size,
% spacing(2) is dim 2's, spacing(3) is dim 3's.
%
% DERIVATION (from the separable 3D Gaussian; applyHessian3D is this
% function with spacing [1 1 1], divided by sigma^2):
%   G(p1,p2,p3) = 1/((2pi)^1.5 s1 s2 s3) * exp(-(p1^2/2s1^2 + p2^2/2s2^2 + p3^2/2s3^2))
%   d2G/dp1^2   = (p1^2/s1^4 - 1/s1^2) * G
%   d2G/dp1dp2  = (p1*p2/(s1^2*s2^2)) * G
%   (s1,s2,s3 = sigmaPhysical(i)./spacing, i.e. PIXEL sigma per axis, chosen
%   so the PHYSICAL smoothing scale along axis i is sigmaPhysical(i))
% then Dij_phys = Dij_pixel / (spacing_i * spacing_j) (chain rule: physical
% position = pixel position * spacing, so d2/dXphys2 = d2/dXpix2 / spacing^2),
% and finally Lindeberg scale normalisation -- Dii by sigmaPhysical(i)^2,
% Dij (i~=j) by sigmaPhysical(i)*sigmaPhysical(j) -- the standard
% anisotropic-scale-space generalisation, per-entry rather than one shared
% scalar. Reduces exactly to the original single-scalar formula when
% sigma1=sigma2=sigma3.
%
% INPUTS
%   I             : 3D volume, native (anisotropic) voxel grid
%   sigmaPhysical : Gaussian scale, in the SAME physical units as spacing
%                   (e.g. microns). Scalar (applied to all 3 axes) or a
%                   3-element [sigma1 sigma2 sigma3] vector (independent
%                   per axis, same dimension order as spacing/size(I)).
%   spacing       : [spacingDim1 spacingDim2 spacingDim3] -- physical size
%                   of one voxel along each array dimension (i.e.
%                   size(I) order). [1 1 1] gives plain isotropic
%                   pixel-space behaviour, identical to
%                   sigma^2 * applyHessian3D (i.e. hessianEigen3D).
%
% OUTPUTS
%   Dxx,Dxy,Dxz,Dyy,Dyz,Dzz : the six unique 2nd derivatives, in PHYSICAL
%                             curvature units (Lindeberg-normalised)
%
% See also: applyHessian3D, hessianEigen3DAniso, eig3volume

if isscalar(sigmaPhysical)
    sigmaPhysical = sigmaPhysical * [1 1 1];
end

s1 = sigmaPhysical(1) / spacing(1);
s2 = sigmaPhysical(2) / spacing(2);
s3 = sigmaPhysical(3) / spacing(3);

r1 = max(1, round(3*s1));
r2 = max(1, round(3*s2));
r3 = max(1, round(3*s3));

% SEPARABLE evaluation (2026-10-03). Every kernel above is a product of 1D
% factors -- e.g. d2G/dp1^2 = g1''(p1)*g2(p2)*g3(p3), d2G/dp1dp2 =
% g1'(p1)*g2'(p2)*g3(p3) -- and replicate padding commutes with separable
% filtering, so 15 1D passes (3 along dim 3, 6 along dim 2, 6 along dim 1,
% sharing the intermediates) give the same result as the six dense 3D
% convolutions this function used to run, at a fraction of the cost
% (kernel taps per voxel drop from ~6*(2r+1)^3 to ~15*(2r+1)). Verified
% against the dense kernels in TestHessian3DAniso. The 1D factors are
% moment-matched rather than the bare sampled formulas -- see localGauss1D.
[g1,h1,k1] = localGauss1D(s1, r1);
[g2,h2,k2] = localGauss1D(s2, r2);
[g3,h3,k3] = localGauss1D(s3, r3);

A0 = localConv(I, g3, 3);  A1 = localConv(I, h3, 3);  A2 = localConv(I, k3, 3);

B00 = localConv(A0, g2, 2);  B01 = localConv(A0, h2, 2);  B02 = localConv(A0, k2, 2);
B10 = localConv(A1, g2, 2);  B11 = localConv(A1, h2, 2);
B20 = localConv(A2, g2, 2);
clear A0 A1 A2

D11_pix = localConv(B00, k1, 1);
D12_pix = localConv(B01, h1, 1);
D22_pix = localConv(B02, g1, 1);
D13_pix = localConv(B10, h1, 1);
D23_pix = localConv(B11, g1, 1);
D33_pix = localConv(B20, g1, 1);
clear B00 B01 B02 B10 B11 B20

% Physical-unit rescaling (chain rule) + PER-ENTRY Lindeberg scale
% normalisation: Dii by sigmaPhysical(i)^2, Dij by sigmaPhysical(i)*sigmaPhysical(j).
Dxx = sigmaPhysical(1)^2 * D11_pix / (spacing(1)*spacing(1));
Dyy = sigmaPhysical(2)^2 * D22_pix / (spacing(2)*spacing(2));
Dzz = sigmaPhysical(3)^2 * D33_pix / (spacing(3)*spacing(3));
Dxy = sigmaPhysical(1)*sigmaPhysical(2) * D12_pix / (spacing(1)*spacing(2));
Dxz = sigmaPhysical(1)*sigmaPhysical(3) * D13_pix / (spacing(1)*spacing(3));
Dyz = sigmaPhysical(2)*sigmaPhysical(3) * D23_pix / (spacing(2)*spacing(3));
end

% =========================================================================
function [g, h, k] = localGauss1D(s, r)
% 1D factors of the 3D kernels, with w = exp(-p^2/2s^2) on p = -r..r:
%   g = a*w            smoothing factor        sum(g)     = 1
%   h = b*p.*w         1st-derivative factor   sum(p.*h)  = 1  (conv gives -f';
%                      the sign cancels in every product used)
%   k = (c + d*p.^2).*w  2nd-derivative factor sum(k) = 0, sum(p.^2.*k) = 2
% MOMENT-MATCHED (2026-10-03). These are the analytic shapes (the continuum
% has a = 1/(sqrt(2pi)s), b = a/s^2, c = -a/s^2, d = a/s^4), but with the
% coefficients fitted so the DISCRETE moments are exact: the Hessian is
% then exact on any cubic, whatever s and r. The analytic coefficients,
% sampled and truncated at r = round(3s), gave sum(k) ~= 0 -- a DC leak:
% Dxx of a uniform region was -0.4% (s=1) to -2% (s=4) of the response to
% unit curvature, -10% at s=0.5, so the in-plane L2 of a bright sheet read
% negative (L2/L3 = 0.01-0.02 isotropic, 0.2-0.4 with a sub-pixel Z sigma
% on a 3x anisotropic grid) and biased the polarity gate and Ra. They also
% gave sum(p.^2.*k)/2 ~= sum(p.*h)^2 (0.89 vs 0.96 at s=4, 1.40 vs 0.77 at
% s=0.5), so diagonal and cross derivatives had different gains and the
% eigenvalues depended on orientation. Enlarging r to 4s would cut the
% truncation part (~20x for s>=2) but not the sampling error at sub-pixel
% s (the usual Z case), and costs 1/3 more taps; this costs nothing.
% For r=1 it reduces to the central differences [1 -2 1] and [1 0 -1]/2.
p = -r:r;
w = exp(-p.^2/(2*s^2));
g = w / sum(w);
h = p .* w;
h = h / sum(p .* h);
M = [sum(w), sum(p.^2.*w); sum(p.^2.*w), sum(p.^4.*w)];
cd = M \ [0; 2];
k = (cd(1) + cd(2)*p.^2) .* w;
end

function B = localConv(A, kern, dim)
sz = ones(1, 3);
sz(dim) = numel(kern);
B = imfilter(A, reshape(kern, sz), 'conv', 'replicate');
end
