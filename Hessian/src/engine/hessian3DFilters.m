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
%   Two ways to handle anisotropic voxels, selected by 'Spacing':
%     - 'Spacing' omitted or [1 1 1] (default): I is treated as an
%       isotropic pixel grid and 'Sigmas' are plain pixel sigmas -- the
%       original behaviour, unchanged, via hessianEigen3D.
%     - 'Spacing' given as [s1 s2 s3] (physical voxel size per array
%       dimension, e.g. microns): I is used on its NATIVE grid, no
%       resampling -- 'Sigmas' are then PHYSICAL sigmas (same units as
%       Spacing), and hessianEigen3DAniso does the per-axis pixel-sigma +
%       physical-unit rescaling (see its header). This avoids the cost of
%       isotropic resampling -- on a real 76-slice, 5.3x anisotropic ER
%       volume, resampling inflated Z to 404 slices and cost ~450s per
%       multiscale vesselness pass; the native-grid path works on the
%       original 76 slices directly.
%
% SUPPORTED FILTERS (FilterType)
%   'vesselness'   - Frangi vesselness (tubular structures)
%   'plate'        - Genuine 3D sheetness (planar/membrane structures,
%                    e.g. ER cisternae -- has no adequate 2D analogue)
%
% INPUTS
%   I              - 3D volume (numeric), native or isotropic grid (see Spacing)
%
% NAME-VALUE PAIRS
%   'FilterType'   - filter to apply (default: 'vesselness')
%   'Sigmas'       - Gaussian scales, one of two forms:
%                      * numeric vector, e.g. [1 2 3] -- N scalar scale
%                        steps, each broadcast to all 3 axes (pixel units
%                        if Spacing is [1 1 1] (default), else physical
%                        units matching Spacing). Unchanged original form.
%                      * cell array of 3-element vectors, e.g.
%                        {[s1 s2 s3], [s1b s2b s3b]} -- N scale steps, each
%                        with an INDEPENDENT physical scale per axis. Only
%                        meaningful with a non-isotropic Spacing. Added for
%                        severely anisotropic data where one shared scalar
%                        sigma cannot serve both axes well: on real 5.3x
%                        anisotropic ER data, forcing enough Z-sigma for
%                        one non-degenerate Z-pixel (sigma>=dz) forced the
%                        SAME scalar through XY, over-smoothing past the
%                        true tubule width and shifting detection from fine
%                        tubules to coarse blobs. Deliberately a cell array
%                        (never a plain Nx3 numeric matrix) so a 3-element
%                        numeric vector is never ambiguous between "3 scalar
%                        steps" and "one per-axis entry".
%   'Spacing'      - [s1 s2 s3], physical voxel size per array dimension;
%                    default [1 1 1] (isotropic pixel-space behaviour,
%                    unchanged from the original engine)
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
    options.Spacing (1,3) double = [1 1 1]
    options.WhiteOnDark (1,1) logical = true
    options.Precision (1,:) char {mustBeMember(options.Precision,{'single','double'})} = 'single'
    options.Parameters = struct()
end

isAniso = ~isequal(options.Spacing, [1 1 1]);

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
    if iscell(sigmas)
        sigma = sigmas{k};
    else
        sigma = sigmas(k);
    end

    if isAniso
        [L1, L2, L3] = hessianEigen3DAniso(I, sigma, options.Spacing, options.Precision);
    else
        [L1, L2, L3] = hessianEigen3D(I, sigma, options.Precision);
    end

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
