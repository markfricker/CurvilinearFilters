function [L1, L2, L3] = hessianEigen3D(I, sigma, Precision)
%HESSIANEIGEN3D  Hessian eigenvalues (3D)
%
%   [L1, L2, L3] = hessianEigen3D(I, sigma, Precision)
%
%   |L1| <= |L2| <= |L3|
%
% 3D counterpart of hessianEigen2D.m. No eigenvector output -- see
% eig3volume.m for why (orientation is unused downstream, same as the 2D
% engine since 2026-09-27).
%
% INPUTS
%   I         : 3D volume, isotropic voxel grid (resample anisotropic data
%               to isotropic first -- this function has no voxel-spacing
%               input, exactly like regionprops3, and would otherwise
%               distort the eigenvalues the same way flagged for
%               trackMitochondria3d).
%   sigma     : Gaussian scale, in voxels
%   Precision : 'single' (default) or 'double'
%
% See also: hessianEigen2D, applyHessian3D, eig3volume

if nargin < 3
    Precision = 'single';
end

if strcmpi(Precision,'single')
    I = single(I);
else
    I = double(I);
end

[Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3D(I, sigma);

% Scale normalisation (Lindeberg) -- matches hessianEigen2D's sigma^2 scaling
Dxx = sigma^2 * Dxx; Dxy = sigma^2 * Dxy; Dxz = sigma^2 * Dxz;
Dyy = sigma^2 * Dyy; Dyz = sigma^2 * Dyz; Dzz = sigma^2 * Dzz;

[L1, L2, L3] = eig3volume(Dxx, Dxy, Dxz, Dyy, Dyz, Dzz);
end
