function imfMasked = localFrobeniusMask3D(imf, im, sigmas, frobDivision)
%LOCALFROBENIUSMASK3D  Zero out voxels with low multi-scale Hessian energy.
%
%   imfMasked = localFrobeniusMask3D(imf, im, sigmas, frobDivision)
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
%                 from, same size as imf. Must already be on an isotropic
%                 voxel grid -- applyHessian3D has no voxel-spacing input,
%                 same caveat as the rest of this engine.
%   sigmas      - vector of Gaussian scales (voxels) to take the per-voxel
%                 max Frobenius norm over -- pass the SAME sigmas used to
%                 produce imf, so the gate reflects the same scale range
%   frobDivision - bias divisor (default 2, matching Nellie's own default
%                 and the 2D port): frobDivision=2 is a deliberately
%                 permissive gate that only clears the clearly-empty
%                 background, leaving the continuous Rb/beta/structureness
%                 terms inside the response formulas to do the fine
%                 discrimination. Larger frobDivision -> lower threshold
%                 -> more permissive (more voxels survive).
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

frobMax = zeros(size(im), 'single');
for s = sigmas
    [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = applyHessian3D(im, s);
    Dxx = s^2*Dxx; Dxy = s^2*Dxy; Dxz = s^2*Dxz;
    Dyy = s^2*Dyy; Dyz = s^2*Dyz; Dzz = s^2*Dzz;

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
