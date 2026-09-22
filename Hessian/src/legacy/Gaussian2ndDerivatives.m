function [Gxx,Gxy,Gyy] = Gaussian2ndDerivatives(s)
%GAUSSIAN2NDDERIVATIVES  2D Gaussian second-derivative convolution kernels.
%
%   [Gxx,Gxy,Gyy] = Gaussian2ndDerivatives(s)
%
% Dependency of HessianMatrix2D -- was missing from this repo; restored
% 2026-09 from C:\Users\dops0035\Documents\Research\Matlab Working\
% granulo_enhancement\code\utilities\Gaussian2ndDerivatives.m.
%
% See also: HessianMatrix2D

%% Grid coordinates
[x,y]   = ndgrid(-round(3*s):round(3*s));
%% Gaussian 2nd derivatives
Gxx = 1/(2*pi*s^4) * (x.^2/s^2 - 1) .* exp(-(x.^2 + y.^2)/(2*s^2));
Gxy = 1/(2*pi*s^6) * (x .* y)       .* exp(-(x.^2 + y.^2)/(2*s^2));
Gyy = Gxx';
%% End
end
