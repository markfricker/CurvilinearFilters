function [lineMap, dirMap] = steerGaussEnhance(im, sigmas, thetas, precision, mode)
% steerGaussEnhance  Multi-scale steerable 2nd-derivative Gaussian ridge filter.
%
% Freeman & Adelson (1991) steerable 2nd-derivative-of-Gaussian filter, as
% used in:
%   Sten et al. (2024) "A Ridge-Based Detection Algorithm with Filament
%   Overlap Identification for 2D Mycelium Network Analysis."
%   DOI: 10.1016/j.ecoinf.2024.102670
%
% The response is the maximum over all scales and orientations; dirMap holds
% the orientation (degrees) that produced the maximum at each pixel.
%
% This is an alternative to AGlineDetectorSteerableConv2 (SOAGK): it uses a
% simpler isotropic Gaussian (no anisotropy factor rho) and is appropriate
% when filament cross-sections are approximately circular.
%
% Inputs
%   im        - 2-D single/double image, normalised [0,1]
%   sigmas    - vector of scales (pixels); default [2 3 4]
%   thetas    - vector of orientations (degrees); default 0:15:345.
%               Ignored in 'analytic' mode.
%   precision - 'single' (default) or 'double' output class
%   mode      - 'sampled' (default): maximum over the listed thetas. Output
%                   matches the original steerGaussFilterOrder2-based
%                   implementation (MycNetAnalysis) to floating-point
%                   rounding.
%               'analytic': exact maximum over ALL orientations, i.e. the
%                   principal eigenvalue of the (scaled) Hessian at each
%                   scale. Slightly larger responses than 'sampled' where
%                   the true orientation falls between sampled angles;
%                   dirMap is continuous, in [0,180).
%
% Outputs
%   lineMap - max ridge response over all scales/orientations, [0,1]
%   dirMap  - orientation (degrees) of max response, same size as im. The
%             direction of maximum curvature, i.e. ACROSS the ridge (the
%             feature normal), +ve anticlockwise.
%
% IMPLEMENTATION: steerability means the oriented response at any theta is
%   J(theta) = cos^2(t)*Ixx + sin^2(t)*Iyy - 2*cos(t)*sin(t)*Ixy,  t = -theta
% so the three basis responses Ixx/Ixy/Iyy are computed ONCE per scale (by
% separable 1-D convolutions) and each orientation is a cheap weighted sum.
% The original implementation recomputed three full 2-D convolutions per
% (scale, orientation) pair. Kernels, normalisation (g0 = exp(-r^2/2s^2) /
% (s*sqrt(2*pi)), no scale normalisation), support (+/- floor(4*sigma)) and
% 'replicate' padding are kept identical so 'sampled' output is unchanged.
%
% POLARITY (fixed 2026-09): the 2nd-derivative kernel is negative at its
% centre, so the raw response is NEGATIVE at bright ridges. It is negated
% before the max/clamp so bright-on-dark ridges give positive output,
% matching the convention used elsewhere in this codebase (e.g.
% FrangiFilter2D's BlackWhite=false). Verified against real bright-on-dark
% ER tubule data.
%
% See also: AGlineDetectorSteerableConv2, applyHessian2D, eig2image

if nargin < 2 || isempty(sigmas),    sigmas    = [2 3 4];      end
if nargin < 3 || isempty(thetas),    thetas    = 0:15:345;     end
if nargin < 4 || isempty(precision), precision = 'single';     end
if nargin < 5 || isempty(mode),      mode      = 'sampled';    end

castfun = str2func(precision);
I = double(im);   % basis responses always computed in double, as before

lineMap = -inf(size(I), precision);
dirMap  = zeros(size(I), precision);

switch lower(mode)
    case 'sampled'
        % Steering weights: J(:,theta) = [Ixx Iyy Ixy] * W(:,theta), all
        % orientations in one matrix product per pixel chunk. max(...,[],2)
        % returns the FIRST maximising theta and the cross-scale update uses
        % strict '>', so ties resolve to the first (scale, theta) exactly as
        % max(respAll,[],3) did in the original stacked version.
        t = -thetas(:).' * (pi/180);
        W = [cos(t).^2; sin(t).^2; -2*cos(t).*sin(t)];
        thetasC = castfun(thetas(:));
        nPix  = numel(I);
        chunk = 2^16;   % pixels per block -- bounds memory to chunk*nTheta doubles
        for is = 1:numel(sigmas)
            [Ixx, Ixy, Iyy] = localBasis(I, sigmas(is));
            for i0 = 1:chunk:nPix
                r = i0 : min(i0+chunk-1, nPix);
                % Negated -- see POLARITY note above.
                J = -castfun([Ixx(r).' Iyy(r).' Ixy(r).'] * W);
                [mx, ix] = max(J, [], 2);
                m = mx > lineMap(r).';
                rm = r(m);
                lineMap(rm) = mx(m);
                dirMap(rm)  = thetasC(ix(m));
            end
        end

    case 'analytic'
        % J(t) = v'*M*v with v = [cos t; sin t], M = [Ixx -Ixy; -Ixy Iyy].
        % max_t(-J) = -lambda_min(M), attained along the lambda_min
        % eigenvector.
        for is = 1:numel(sigmas)
            [Ixx, Ixy, Iyy] = localBasis(I, sigmas(is));
            half = 0.5*(Ixx - Iyy);
            R    = castfun(-(0.5*(Ixx + Iyy) - sqrt(half.^2 + Ixy.^2)));
            m = R > lineMap;
            lineMap(m) = R(m);
            % lambda_max axis angle: 0.5*atan2(2*(-Ixy), Ixx-Iyy); lambda_min
            % is perpendicular. theta = -t, reported in [0,180).
            tMin = 0.5*atan2(-2*Ixy, Ixx - Iyy) + pi/2;
            th   = castfun(mod(-tMin*(180/pi), 180));
            dirMap(m) = th(m);
        end

    otherwise
        error('steerGaussEnhance:mode', ...
            'Unknown mode ''%s'' (use ''sampled'' or ''analytic'').', mode);
end

% Clamp negatives -- after the sign flip, negative means the response
% favours a dark ridge, i.e. no bright ridge at that pixel.
lineMap(lineMap < 0) = castfun(0);

% Normalise to [0,1]
mx = max(lineMap(:));
if mx > 0
    lineMap = lineMap ./ mx;
end

end

% -------------------------------------------------------------------------
function [Ixx, Ixy, Iyy] = localBasis(I, sigma)
% Separable 2nd-derivative-of-Gaussian basis responses. Equivalent to
% correlating with the 2-D kernels
%   G2a = g0.*(x.^2/s^4 - 1/s^2),  G2b = g0.*x.*y/s^4,  G2c = G2a'
%   g0  = exp(-(x.^2+y.^2)/(2 s^2)) / (s*sqrt(2*pi))
% with x along columns, 'replicate' padding (replicate padding is itself
% separable, so the two-pass result equals the 2-D one up to rounding).
W = max(1, floor(4*sigma));
x = -W:W;
g   = exp(-x.^2 / (2*sigma^2));               % un-normalised 1-D Gaussian
gN  = g / (sigma*sqrt(2*pi));                 % carries g0's normalisation
d2  = gN .* (x.^2/sigma^4 - 1/sigma^2);       % normalised 2nd-derivative
d1N = gN .* x / sigma^2;                      % normalised 1st-derivative
d1  = g  .* x / sigma^2;                      % un-normalised 1st-derivative

Ixx = imfilter(imfilter(I, d2,  'replicate'), g.',  'replicate');
Iyy = imfilter(imfilter(I, g,   'replicate'), d2.', 'replicate');
Ixy = imfilter(imfilter(I, d1N, 'replicate'), d1.', 'replicate');
end
