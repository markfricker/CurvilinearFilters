function [response, scale, direction] = hessian2DFilters(I, options)
%HESSIAN2DFILTERS Multiscale Hessian-based filtering (2D)
%
%   [R, SCALE, DIR] = hessian2DFilters(I, Name, Value, ...)
%
% DESCRIPTION
%   hessian2DFilters implements a clean, stateless, multiscale
%   Hessian-based filtering framework for 2D images.
%
%   The algorithm:
%     1) builds a Gaussian scale space over the specified SIGMAS;
%     2) computes the Hessian matrix at each scale;
%     3) extracts Hessian eigenvalues and eigenvectors;
%     4) evaluates a scalar, per-scale response function;
%     5) aggregates responses by MAX-over-scales.
%
%   This engine supports only Hessian-compatible filters
%
% SUPPORTED FILTERS (FilterType)
%   'vesselness'   - Frangi vesselness (tubular structures)
%   'ridge'        - Sato-style ridge detector (line-likeness)
%   'blob'         - Blobness (determinant-of-Hessian-like)
%   'plate'        - Plateness (2D analogue for thick elongated regions)
%
% INPUTS
%   I              - 2D grayscale image (numeric)
%
% NAME–VALUE PAIRS
%   'FilterType'   - filter to apply (default: 'vesselness')
%   'Sigmas'       - vector of Gaussian scales (default: auto)
%   'WhiteOnDark'  - true for bright structures on dark background (default: true)
%   'Precision'    - 'single' or 'double' (default: 'single')
%   'Parameters'   - struct of filter-specific parameters
%
% OUTPUTS
%   R              - response image (same size as I)
%   SCALE          - index of scale where maximum response occurred
%   DIR            - local orientation (radians), from minor eigenvector
%
% DESIGN CONTRACT
%   - max-over-scales aggregation only
%   - stateless per scale
%

arguments
    I (:,:) {mustBeNumeric}
    options.FilterType (1,:) char = 'vesselness'
    options.Sigmas = []
    options.WhiteOnDark (1,1) logical = true
    options.Precision (1,:) char {mustBeMember(options.Precision,{'single','double'})} = 'single'
    options.Parameters = struct()
end

% Presets
options = hessian2DPresets(I, options);
sigmas = options.Sigmas;

% Outputs
response  = zeros(size(I), 'like', I);
scale     = zeros(size(I), 'uint16');
direction = zeros(size(I), 'like', I);

% Response function
switch lower(options.FilterType)

    case 'vesselness'
        responseFcn = @(L1,L2) vesselnessResponse( ...
            L1, L2, options.Parameters.beta, options.Parameters.c);

    case 'ridge'
        responseFcn = @(L1,L2) ridgeResponse( ...
            L1, L2, options.Parameters.alpha);

    case 'blob'
        responseFcn = @(L1,L2) blobnessResponse(L1, L2);

    case 'plate'
        responseFcn = @(L1,L2) platenessResponse( ...
            L1, L2, options.Parameters.alpha);

    otherwise
        error('FilterType "%s" not supported in hessian2DFilters.', ...
              options.FilterType);
end

% Multiscale loop
for k = 1:numel(sigmas)
    sigma = sigmas(k);

    [L1, L2, Vx, Vy] = hessianEigen2D(I, sigma, options.Precision);

    R = responseFcn(L1, L2);

    % Polarity (standard Frangi/Sato semantics)
    if options.WhiteOnDark
        R(L2 >= 0) = 0;
    else
        R(L2 <= 0) = 0;
    end

    % Max-over-scales
    mask = R > response;
    response(mask) = R(mask);
    scale(mask) = k;

    % Tangent direction (perpendicular to Hessian normal)
    direction(mask) = atan2(-Vy(mask), Vx(mask));
end
end
