function [imf,Vx,Vy] = NeuritenessFilter2D(im,sigma)
%%  NeuriteneesFilter - computing the neuritenees
%   
%   REFERENCE:
%       E. Meijering et al. Design and Validation of a Tool for Neurite 
%       Tracing and Analysis in Fluorescence Microscopy Images
%       Cytometry Part A, 58A, 167–176, 2004
%
%   INPUT:
%       im      - gray image
%       sigma   - sigma factor
%
%   OUTPUT:
%       imf     - filtered image
%
%   USAGE:
%
%   AUTHOR:
%       Boguslaw Obara, http://boguslawobara.net/
%
%   VERSION:
%       0.1 - 25/04/2010 First implementation

%% Cost Map
[L1,L2,V1,V2,V3,V4] = HessianVectorField(im,sigma);
% Modified Hessian
L1t = L1;
L2t = L2;
alfa = -1/3;
L1 = L1t + alfa*L2t;
L2 = L2t + alfa*L1t;
% Sort L1s > L2s
index = abs(L1)<abs(L2);
L1s = L1; 
L2s = L2;
L1s(index) = L2(index);
L2s(index) = L1(index);
%V1s = V1; 

V3s = V3;
%V2s = V2; 
V4s = V4;
%V1s(index) = V3(index);
V3s(index) = V1(index);
%V2s(index) = V4(index);
V4s(index) = V2(index);

% Neuriteness
Lmin = min(L1s(:));
imf = zeros(size(L1));
imf(L1s<0) = L1s(L1s<0)./Lmin;
% Eigenvectors
%Vy = V4s; Vx = V3s;
Vy = V4; Vx = V3;

% name = strcat(name,'Neurite','.png');
% saveFunction(imf,'2D',name);
end