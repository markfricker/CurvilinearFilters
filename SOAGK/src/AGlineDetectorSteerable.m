function [lineMap, dirMap] = AGlineDetectorSteerable(im, sigmas, thetas, rhos, precision)
% AGlineDetectorSteerable
%
% Optimized anisotropic Gaussian line detector using steerable SOAGK
%
% Inputs
%   im        : grayscale image in [0,1]
%   sigmas    : vector of scales
%   thetas    : vector of orientations [0,pi)
%   rhos      : anisotropy factors >= 1
%   precision : 'single' or 'double'
%
% Outputs
%   lineMap : max line response
%   dirMap  : orientation map (radians)

% ----------------------------------------------------
% Precision handling
% ----------------------------------------------------
if nargin < 5
    precision = 'single';
end

castfun = str2func(precision);
im      = castfun(im);

lineMap = zeros(size(im), precision);
dirMap  = zeros(size(im), precision);

% ----------------------------------------------------
% Main loop
% ----------------------------------------------------
for s = 1:numel(sigmas)
    sigma = sigmas(s);
    sigma2 = sigma^2;

    for r = 1:numel(rhos)
        rho = rhos(r);

        % ------------------------------------------------
        % Create steerable basis filters (once!)
        % ------------------------------------------------
        params.rho   = rho;
        params.theta = 0; % axis-aligned
        params.size  = round(9*sigma*max(1,1/rho));
        if mod(params.size,2)==0
            params.size = params.size + 1;
        end

        [Gxx, Gyy, Gxy] = createSOAGKSteerableBasis(sigma, params, precision);

        % ------------------------------------------------
        % Convolve once per basis
        % ------------------------------------------------
        Rxx = imfilter(im, Gxx, 'replicate');
        Ryy = imfilter(im, Gyy, 'replicate');
        Rxy = imfilter(im, Gxy, 'replicate');

        % ------------------------------------------------
        % Steer responses
        % ------------------------------------------------
        for t = 1:numel(thetas)
            th = thetas(t);
            c  = cos(th);
            s  = sin(th);

            response = ...
                (c*c) * Rxx + ...
                (s*s) * Ryy + ...
                (2*c*s) * Rxy;

            % Normalize (zero-mean line response)
            response = -response;

            mask = response > lineMap;
            lineMap(mask) = response(mask);
            dirMap(mask)  = th;
        end
    end
end

% Normalize lineMap to [0,1]
lineMap = lineMap - min(lineMap(:));
lineMap = lineMap / max(lineMap(:));
end
