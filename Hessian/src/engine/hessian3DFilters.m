function [response, scale] = hessian3DFilters(I, options)
%HESSIAN3DFILTERS Multiscale Hessian-based filtering (3D)
%
%   [R, SCALE] = hessian3DFilters(I, Name, Value, ...)
%
% DESCRIPTION
%   3D counterpart of hessian2DFilters.m. The algorithm:
%     1) builds a Gaussian scale space over the specified SIGMAS;
%     2) computes the 3D Hessian matrix at each scale;
%     3) extracts the three Hessian eigenvalues (eig3volume, no
%        eigenvectors -- see its header for why);
%     4) evaluates a scalar, per-scale response function using all three
%        eigenvalues (genuinely different structures for vesselness vs
%        plateness -- see those functions' headers);
%     5) aggregates responses by MAX-over-scales.
%
%   I MUST be on an isotropic voxel grid (resample anisotropic data
%   first) -- eig3volume/regionprops3-style eigenvalue-based measures have
%   no voxel-spacing input and are distorted by anisotropic voxels
%   exactly as documented for trackMitometer3d.m.
%
% SUPPORTED FILTERS (FilterType)
%   'vesselness'   - Frangi vesselness (tubular structures)
%   'plate'        - Genuine 3D sheetness (planar/membrane structures,
%                    e.g. ER cisternae -- has no adequate 2D analogue)
%
% INPUTS
%   I              - 3D volume (numeric), isotropic voxel grid
%
% NAME-VALUE PAIRS
%   'FilterType'   - filter to apply (default: 'vesselness')
%   'Sigmas'       - vector of Gaussian scales, in voxels (default: auto)
%   'WhiteOnDark'  - true for bright structures on dark background (default: true)
%   'Precision'    - 'single' or 'double' (default: 'single')
%   'Parameters'   - struct of filter-specific parameters (alpha, beta, c)
%
% OUTPUTS
%   R              - response volume (same size as I)
%   SCALE          - index of scale where maximum response occurred
%
% DESIGN CONTRACT
%   - max-over-scales aggregation only, stateless per scale
%   - no orientation output (mirrors the 2D engine since 2026-09-27 --
%     nothing downstream reads it; 3D eigenvectors are also the genuinely
%     hard part of this problem, deliberately not built until something
%     needs them)
%
% See also: hessian2DFilters, hessianEigen3D, vesselnessResponse3D,
%   platenessResponse3D

arguments
    I (:,:,:) {mustBeNumeric}
    options.FilterType (1,:) char = 'vesselness'
    options.Sigmas = []
    options.WhiteOnDark (1,1) logical = true
    options.Precision (1,:) char {mustBeMember(options.Precision,{'single','double'})} = 'single'
    options.Parameters = struct()
end

options = hessian3DPresets(I, options);
sigmas = options.Sigmas;

response = zeros(size(I), 'like', I);
scale    = zeros(size(I), 'uint16');

switch lower(options.FilterType)
    case 'vesselness'
        responseFcn = @(L1,L2,L3) vesselnessResponse3D( ...
            L1, L2, L3, options.Parameters.alpha, options.Parameters.beta, options.Parameters.c);
    case 'plate'
        responseFcn = @(L1,L2,L3) platenessResponse3D( ...
            L1, L2, L3, options.Parameters.alpha, options.Parameters.beta, options.Parameters.c);
    otherwise
        error('FilterType "%s" not supported in hessian3DFilters.', options.FilterType);
end

for k = 1:numel(sigmas)
    sigma = sigmas(k);

    [L1, L2, L3] = hessianEigen3D(I, sigma, options.Precision);

    R = responseFcn(L1, L2, L3);

    % Polarity gate (standard Frangi/Sato semantics, mirrors
    % hessian2DFilters.m exactly): the response functions are
    % polarity-symmetric, so one eigenvalue pass serves either sign.
    if options.WhiteOnDark
        R(L2 >= 0 | L3 >= 0) = 0;
    else
        R(L2 <= 0 | L3 <= 0) = 0;
    end

    mask = R > response;
    response(mask) = R(mask);
    scale(mask) = k;
end
end
