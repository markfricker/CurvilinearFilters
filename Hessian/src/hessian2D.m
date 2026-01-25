function [response, scale, direction] = hessian2D(I, sigmas, responseFcn)
% hessian2D  Clean multiscale Hessian-based filter (2D)
%
% Inputs:
%   I           - 2D image (single or double)
%   sigmas      - vector of scales
%   responseFcn - handle @(l1,l2) -> scalar response
%
% Outputs:
%   response  - max response over scales
%   scale     - scale index of max response
%   direction - orientation (minor eigenvector)

I = single(I);
response = zeros(size(I),'single');
scale = zeros(size(I),'uint16');
direction = zeros(size(I),'single');

for k = 1:numel(sigmas)
    sigma = sigmas(k);

    % Hessian
    [Dxx,Dxy,Dyy] = Hessian2D(I, sigma);
    Dxx = sigma^2 * Dxx;
    Dxy = sigma^2 * Dxy;
    Dyy = sigma^2 * Dyy;

    % Eigenvalues / vectors
    [l1,l2,vx,vy] = eig2image(Dxx,Dxy,Dyy);

    % Response at this scale
    R = responseFcn(l1,l2);

    % Polarity: bright-on-dark
    R(l2 >= 0) = 0;

    % Max-over-scales
    mask = R > response;
    response(mask) = R(mask);
    scale(mask) = k;
    direction(mask) = atan2(vx(mask), vy(mask));
end
end
