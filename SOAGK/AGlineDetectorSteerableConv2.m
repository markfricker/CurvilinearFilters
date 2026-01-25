function [lineMap, dirMap] = AGlineDetectorSteerableConv2( ...
        im, sigmas, thetas, rhos, precision)
% AGlineDetectorSteerableConv2
%
% Fast anisotropic Gaussian line detector using steerable SOAGK
% and conv2 (datatype-safe, low memory).
%
% Inputs
%   im        : grayscale image [0,1]
%   sigmas    : vector of scales
%   thetas    : orientations in [0, pi)
%   rhos      : anisotropy factors >= 1
%   precision : 'single' or 'double'
%
% Outputs
%   lineMap : line strength map
%   dirMap  : orientation map (radians)

% ----------------------------
% Precision handling
% ----------------------------
if nargin < 5
    precision = 'single';
end
castfun = str2func(precision);

im = castfun(im);

lineMap = zeros(size(im), precision);
dirMap  = zeros(size(im), precision);

% ----------------------------
% Main loops
% ----------------------------
for is = 1:numel(sigmas)
    sigma = sigmas(is);

    for ir = 1:numel(rhos)
        rho = rhos(ir);

        % --- build steerable basis once ---
        params.rho = rho;
        params.size = round(9 * sigma * max(1, 1/rho));
        if mod(params.size,2)==0
            params.size = params.size + 1;
        end

        % [Gxx, Gyy, Gxy] = createSOAGKSteerableBasisConv2( ...
        %                     sigma, params, precision);

        [Gxx, Gyy, Gxy] = getSOAGKSteerableBasisCached( ...
                    sigma, rho, precision);

        % --- pad to avoid boundary issues ---
        pad = floor(size(Gxx,1)/2);
        imp = padarray(im,[pad pad],'replicate','both');

        % --- convolve once per basis ---
        Rxx = conv2(im, Gxx, 'same');
        Ryy = conv2(im, Gyy, 'same');
        Rxy = conv2(im, Gxy, 'same');

        % --- steer responses ---
        for it = 1:numel(thetas)
            th = thetas(it);
            c  = cos(th);
            s  = sin(th);

            response = ...
                (c*c) * Rxx + ...
                (s*s) * Ryy + ...
                (2*c*s) * Rxy;

            % We want line-like (negative second derivative)
            response = -response;

            mask = response > lineMap;
            lineMap(mask) = response(mask);
            dirMap(mask)  = th;
        end
    end
end

% ----------------------------
% Normalize output to [0,1]
% ----------------------------
lineMap = lineMap - min(lineMap(:));
mx = max(lineMap(:));
if mx > 0
    lineMap = lineMap / mx;
end
end
