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
% See also: applyHessian2D, hessianEigen3D, eig3volume

if nargin < 2, Sigma = 1; end

r = round(3*Sigma);
[X,Y,Z] = ndgrid(-r:r);

g = @(X,Y,Z) exp(-(X.^2+Y.^2+Z.^2)/(2*Sigma^2));

DGaussxx = 1/(2*pi*Sigma^4) * (X.^2/Sigma^2 - 1) .* g(X,Y,Z);
DGaussyy = permute(DGaussxx, [2 1 3]);
DGausszz = permute(DGaussxx, [3 2 1]);

DGaussxy = 1/(2*pi*Sigma^6) * (X.*Y) .* g(X,Y,Z);
DGaussxz = 1/(2*pi*Sigma^6) * (X.*Z) .* g(X,Y,Z);
DGaussyz = 1/(2*pi*Sigma^6) * (Y.*Z) .* g(X,Y,Z);

Dxx = imfilter(I,DGaussxx,'conv','replicate');
Dyy = imfilter(I,DGaussyy,'conv','replicate');
Dzz = imfilter(I,DGausszz,'conv','replicate');
Dxy = imfilter(I,DGaussxy,'conv','replicate');
Dxz = imfilter(I,DGaussxz,'conv','replicate');
Dyz = imfilter(I,DGaussyz,'conv','replicate');
end
