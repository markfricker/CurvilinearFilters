function imfMasked = localFrobeniusMask3D(imf, im, sigmas, frobDivision, spacing)
%LOCALFROBENIUSMASK3D  Zero out voxels with low multi-scale Hessian energy.
%
%   imfMasked = localFrobeniusMask3D(imf, im, sigmas, frobDivision, spacing)
%
% 3D counterpart of funcEnhanceDispatcher.m's local function
% localFrobeniusMask (AnalyzERproject_sandbox/core/dispatchers/). Same
% technique Nellie uses to gate its own Frangi filter: threshold the
% Hessian's Frobenius norm and hard-zero the enhanced image wherever it
% falls below that -- vesselnessResponse3D/platenessResponse3D are
% continuous functions that only asymptotically approach zero for
% background, so left alone they leave faint background residue visible
% after normalisation. This reproduces the gate as a separate step,
% exactly mirroring the 2D fix that measurably cleaned up real ER vesselness
% output (see [[nellie-comparison-and-auto-threshold-2026-09-22]]).
%
% INPUTS
%   imf         - the enhanced (vesselness/plateness/etc.) volume to mask,
%                 any size/class
%   im          - the ORIGINAL (pre-enhance) volume the Hessian is computed
%                 from, same size as imf, on the SAME grid (native or
%                 isotropic) that produced imf -- see `spacing`.
%   sigmas      - Gaussian scales to take the per-voxel max Frobenius norm
%                 over -- pass the SAME sigmas (same form, same units:
%                 pixel if spacing=[1 1 1], physical otherwise) used to
%                 produce imf, so the gate reflects the same scale range
%                 and axis-weighting hessian3DFilters actually used. Either
%                 a numeric vector (scalar steps, broadcast to all axes) or
%                 a cell array of 3-element per-axis vectors -- see
%                 hessian3DFilters.m's 'Sigmas' for the two forms.
%   frobDivision - bias divisor (default 2, matching Nellie's own default
%                 and the 2D port): frobDivision=2 is a deliberately
%                 permissive gate that only clears the clearly-empty
%                 background, leaving the continuous Rb/beta/structureness
%                 terms inside the response formulas to do the fine
%                 discrimination. Larger frobDivision -> lower threshold
%                 -> more permissive (more voxels survive).
%   spacing     - [s1 s2 s3], physical voxel size per array dimension;
%                 default [1 1 1] (isotropic pixel-space, unchanged
%                 original behaviour via applyHessian3D). Pass the SAME
%                 Spacing given to hessian3DFilters when imf came from its
%                 anisotropic path -- otherwise the gate would pool raw
%                 pixel-space curvatures across axes with very different
%                 physical meaning (e.g. a coarse Z pixel's curvature
%                 looking artificially small next to a fine XY pixel's),
%                 systematically under-weighting real axial structure.
%
% OUTPUT
%   imfMasked   - imf with background voxels (Frobenius norm below the
%                 auto threshold) hard-zeroed
%
% DEPENDENCIES
%   applyHessian3D (this repo)
%   globalThresholdFast (Segmentation_sandbox/src/) -- already dimension-
%     agnostic (pure histogram/threshold, no spatial assumption), reused
%     here on a 3D Frobenius-norm volume with no changes needed
%
% See also: hessian3DFilters, vesselnessResponse3D, platenessResponse3D

if nargin < 4 || isempty(frobDivision)
    frobDivision = 2;
end
if nargin < 5 || isempty(spacing)
    spacing = [1 1 1];
end
isAniso = ~isequal(spacing, [1 1 1]);

frobMax = zeros(size(im), 'single');
for k = 1:numel(sigmas)
    if iscell(sigmas)
        s = sigmas{k};
    else
        s = sigmas(k);
    end
    if isAniso
        [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3DAniso(im, s, spacing);
    else
        [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3D(im, s);
        Dxx = s^2*Dxx; Dxy = s^2*Dxy; Dxz = s^2*Dxz;
        Dyy = s^2*Dyy; Dyz = s^2*Dyz; Dzz = s^2*Dzz;
    end

    % Frobenius norm of a symmetric matrix: sqrt(sum of squares of all
    % entries) = sqrt(diagonal^2 sum + 2*off-diagonal^2 sum) -- the exact
    % 3D generalisation of the 2D gate's sqrt(Dxx^2+2*Dxy^2+Dyy^2).
    F = sqrt(Dxx.^2 + Dyy.^2 + Dzz.^2 + 2*Dxy.^2 + 2*Dxz.^2 + 2*Dyz.^2);
    frobMax = max(frobMax, single(F));
end

mask = globalThresholdFast(frobMax, 'method', 'triangleOtsu', ...
    'bias', 1 / frobDivision);

imfMasked = imf;
imfMasked(~mask) = 0;
end
