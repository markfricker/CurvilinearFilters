function [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3D(I,Sigma)
%APPLYHESSIAN3D  3D Hessian via 2nd derivatives of a Gaussian.
%
%   [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3D(I,Sigma)
%
% 3D counterpart of applyHessian2D.m (same author's Kroon-style kernel
% construction, generalised to three dimensions). The six unique second
% partial derivatives of a 3D volume, convolved with a Gaussian kernel of
% the given Sigma (isotropic in voxel units -- resample the volume to
% isotropic voxels first if the physical voxel spacing is anisotropic).
%
% INPUTS
%   I     : 3D volume, class preferably double or single
%   Sigma : Gaussian kernel sigma, in voxels (default 1)
%
% OUTPUTS
%   Dxx,Dxy,Dxz,Dyy,Dyz,Dzz : the six unique 2nd derivatives (the Hessian
%                             is symmetric: Dyx=Dxy, Dzx=Dxz, Dzy=Dyz)
%
% KERNELS (2026-10-03): computed by applyHessian3DAniso with unit spacing,
% i.e. separable, moment-matched 1D factors (see its localGauss1D). The
% dense kernels used here before, (X^2/Sigma^2 - 1)*exp(...) truncated at
% round(3*Sigma), did not sum to zero: a uniform region gave Dxx < 0 in
% proportion to its intensity, so a bright sheet's in-plane L2 read
% negative (L2/L3 ~0.01-0.02 instead of 0).
%
% NORMALISATION: the true derivative of the unit-mass 3D Gaussian (no scale
% normalisation; hessianEigen3D applies Lindeberg's sigma^2). Before
% 2026-10-03 the kernel carried the 2D constant 1/(2*pi*Sigma^4), making
% the output sqrt(2*pi)*Sigma times too large: net sigma^3 normalisation
% after hessianEigen3D, so the scale-normalised curvature of a blurred
% sheet rose monotonically with sigma instead of peaking at its matched
% scale, biasing scale selection coarse and making a fixed Frangi c mean
% something different at every sigma. hessian3DFilters now derives c from
% the data (hessian3DFrangiC) unless one is given.

if nargin < 2, Sigma = 1; end

[Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3DAniso(I, Sigma, [1 1 1]);

% applyHessian3DAniso returns Sigma^2 * true derivative (Lindeberg)
f = 1 / Sigma^2;
Dxx = f*Dxx; Dxy = f*Dxy; Dxz = f*Dxz;
Dyy = f*Dyy; Dyz = f*Dyz; Dzz = f*Dzz;
end
