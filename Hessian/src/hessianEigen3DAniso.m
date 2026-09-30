function [L1, L2, L3] = hessianEigen3DAniso(I, sigmaPhysical, spacing, Precision)
%HESSIANEIGEN3DANISO  Hessian eigenvalues on a native anisotropic voxel grid.
%
%   [L1, L2, L3] = hessianEigen3DAniso(I, sigmaPhysical, spacing, Precision)
%
%   |L1| <= |L2| <= |L3|, in PHYSICAL curvature units (comparable across
%   axes and scales despite anisotropic voxels -- see applyHessian3DAniso's
%   header for the derivation). No resampling of I -- works on the native
%   grid directly.
%
% INPUTS
%   I             : 3D volume, native (anisotropic) voxel grid
%   sigmaPhysical : Gaussian scale, physical units (e.g. microns)
%   spacing       : [spacingDim1 spacingDim2 spacingDim3], same units as
%                   sigmaPhysical (size(I) dimension order)
%   Precision     : 'single' (default) or 'double'
%
% See also: hessianEigen3D (isotropic), applyHessian3DAniso, eig3volume

if nargin < 4
    Precision = 'single';
end

if strcmpi(Precision,'single')
    I = single(I);
else
    I = double(I);
end

% Lindeberg + physical-unit scaling both already applied inside
% applyHessian3DAniso -- unlike hessianEigen3D (isotropic), which applies
% Lindeberg's sigma^2 scaling itself after a "bare" applyHessian3D call.
[Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3DAniso(I, sigmaPhysical, spacing);

[L1, L2, L3] = eig3volume(Dxx, Dxy, Dxz, Dyy, Dyz, Dzz);
end
