function [lineMap, dirMap] = steerGaussEnhance(im, sigmas, thetas, precision)
% steerGaussEnhance  Multi-scale steerable 2nd-derivative Gaussian ridge filter.
%
% Uses the Freeman & Adelson (1991) steerable filter formulation implemented
% by steerGaussFilterOrder2 (Pang, Tufts 2013), as used in:
%   Sten et al. (2024) "A Ridge-Based Detection Algorithm with Filament
%   Overlap Identification for 2D Mycelium Network Analysis."
%   DOI: 10.1016/j.ecoinf.2024.102670
%
% The response is the maximum over all scales and orientations; dirMap holds
% the orientation (degrees, [0,360)) that produced the maximum at each pixel.
%
% This is an alternative to AGlineDetectorSteerableConv2 (SOAGK): it uses a
% simpler isotropic Gaussian (no anisotropy factor rho) and is appropriate
% when filament cross-sections are approximately circular.
%
% Inputs
%   im        - 2-D single/double image, normalised [0,1]
%   sigmas    - vector of scales (pixels); default [2 3 4]
%   thetas    - vector of orientations (degrees); default 0:15:345
%   precision - 'single' (default) or 'double'
%
% Outputs
%   lineMap - single, max ridge response over all scales/orientations [0,1]
%   dirMap  - single, orientation (degrees) of max response, same size as im
%
% Dependencies
%   steerGaussFilterOrder2  (MycNetAnalysis/src/fcn/external/ must be on path)
%
% See also: AGlineDetectorSteerableConv2, steerGaussFilterOrder2

if nargin < 2 || isempty(sigmas),    sigmas    = [2 3 4];      end
if nargin < 3 || isempty(thetas),    thetas    = 0:15:345;     end
if nargin < 4 || isempty(precision), precision = 'single';     end

castfun = str2func(precision);
im      = castfun(im);

nS = numel(sigmas);
nT = numel(thetas);

respAll = zeros([size(im), nS*nT], precision);

idx = 0;
for is = 1:nS
    for it = 1:nT
        idx = idx + 1;
        respAll(:,:,idx) = castfun( ...
            steerGaussFilterOrder2(im, thetas(it), sigmas(is), false));
    end
end

% Maximum response and corresponding orientation index
[lineMap, maxIdx] = max(respAll, [], 3);

% Recover orientation from flat index
thetaIdx = mod(maxIdx - 1, nT) + 1;
dirMap   = castfun(thetas(thetaIdx));

% Clamp negatives (2nd-derivative filter responds negatively to dark ridges)
lineMap(lineMap < 0) = castfun(0);

% Normalise to [0,1]
mx = max(lineMap(:));
if mx > 0
    lineMap = lineMap ./ mx;
end

end
