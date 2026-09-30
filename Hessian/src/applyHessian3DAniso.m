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
% DERIVATION (from the separable 3D Gaussian, not a rescaled copy of
% applyHessian3D.m's isotropic kernel -- that file uses its own,
% independently self-consistent normalisation convention that does not
% generalise correctly to anisotropic sigma by simple substitution):
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
%                   size(I) order). Pass [1 1 1] to recover plain
%                   isotropic pixel-space behaviour (a DIFFERENT,
%                   independently self-consistent normalisation from
%                   applyHessian3D.m -- use that function directly for the
%                   isotropic case; this one is for genuinely anisotropic
%                   spacing).
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
[P1,P2,P3] = ndgrid(-r1:r1, -r2:r2, -r3:r3);

g = exp(-(P1.^2/(2*s1^2) + P2.^2/(2*s2^2) + P3.^2/(2*s3^2)));
normConst = 1 / ((2*pi)^1.5 * s1 * s2 * s3);
G = normConst * g;

DGauss11 = (P1.^2/s1^4 - 1/s1^2) .* G;
DGauss22 = (P2.^2/s2^4 - 1/s2^2) .* G;
DGauss33 = (P3.^2/s3^4 - 1/s3^2) .* G;
DGauss12 = (P1.*P2 / (s1^2*s2^2)) .* G;
DGauss13 = (P1.*P3 / (s1^2*s3^2)) .* G;
DGauss23 = (P2.*P3 / (s2^2*s3^2)) .* G;

D11_pix = imfilter(I, DGauss11, 'conv', 'replicate');
D22_pix = imfilter(I, DGauss22, 'conv', 'replicate');
D33_pix = imfilter(I, DGauss33, 'conv', 'replicate');
D12_pix = imfilter(I, DGauss12, 'conv', 'replicate');
D13_pix = imfilter(I, DGauss13, 'conv', 'replicate');
D23_pix = imfilter(I, DGauss23, 'conv', 'replicate');

% Physical-unit rescaling (chain rule) + PER-ENTRY Lindeberg scale
% normalisation: Dii by sigmaPhysical(i)^2, Dij by sigmaPhysical(i)*sigmaPhysical(j).
Dxx = sigmaPhysical(1)^2 * D11_pix / (spacing(1)*spacing(1));
Dyy = sigmaPhysical(2)^2 * D22_pix / (spacing(2)*spacing(2));
Dzz = sigmaPhysical(3)^2 * D33_pix / (spacing(3)*spacing(3));
Dxy = sigmaPhysical(1)*sigmaPhysical(2) * D12_pix / (spacing(1)*spacing(2));
Dxz = sigmaPhysical(1)*sigmaPhysical(3) * D13_pix / (spacing(1)*spacing(3));
Dyz = sigmaPhysical(2)*sigmaPhysical(3) * D23_pix / (spacing(2)*spacing(3));
end
