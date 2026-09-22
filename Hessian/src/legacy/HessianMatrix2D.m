function [Hxx,Hxy,Hyy] = HessianMatrix2D(im,s)
%HESSIANMATRIX2D  2D Hessian components via Gaussian 2nd-derivative kernels.
%
%   [Hxx,Hxy,Hyy] = HessianMatrix2D(im,s)
%
% Dependency of HessianVectorField (called by the legacy-path
% NeuritenessFilter2D) -- was missing from this repo; restored 2026-09
% from C:\Users\dops0035\Documents\Research\Matlab Working\
% granulo_enhancement\code\utilities\HessianMatrix2D.m.
%
% See also: Gaussian2ndDerivatives, HessianVectorField, NeuritenessFilter2D

%% Gaussian 2nd derivatives
[Gxx,Gxy,Gyy] = Gaussian2ndDerivatives(s) ;
%% Hessian Matrix
Hxx = imfilter(im,Gxx,'conv','same','replicate');
Hxy = imfilter(im,Gxy,'conv','same','replicate');
Hyy = imfilter(im,Gyy,'conv','same','replicate');
%% End
end
