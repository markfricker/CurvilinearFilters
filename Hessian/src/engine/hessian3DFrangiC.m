function c = hessian3DFrangiC(I, sigmas, spacing)
%HESSIAN3DFRANGIC  Data-driven Frangi "c" for the 3D Hessian filters.
%
%   c = hessian3DFrangiC(I, sigmas, spacing)
%
% c = max over voxels and scales of the scale-normalised Hessian Frobenius
% norm S = sqrt(L1^2+L2^2+L3^2), divided by 2 -- Frangi et al. (1998)'s
% suggestion of half the maximum Hessian norm, and the same frobDivision=2
% convention localFrobeniusMask3D uses. c sets where the structure term
% (1 - exp(-S^2/2c^2)) saturates, so it depends on the data's intensity and
% derivative scale; a fixed constant (the old default 15) does not carry
% between datasets.
%
% INPUTS
%   I       - 3D volume
%   sigmas  - the SAME 'Sigmas' given to hessian3DFilters (numeric vector or
%             cell array of per-axis 3-vectors)
%   spacing - the SAME 'Spacing' (default [1 1 1])
%
% Costs one Hessian pass per scale (no eigen-decomposition: the Frobenius
% norm of a symmetric matrix is sqrt(sum of squared entries)).
%
% See also: hessian3DFilters, hessian3DPresets, localFrobeniusMask3D

if nargin < 3 || isempty(spacing)
    spacing = [1 1 1];
end
if ~isfloat(I)
    I = single(I);
end

Smax = 0;
for k = 1:numel(sigmas)
    if iscell(sigmas)
        s = sigmas{k};
    else
        s = sigmas(k);
    end
    % hessianEigen3D's convention is sigma^2 * applyHessian3D, which equals
    % applyHessian3DAniso with spacing [1 1 1] -- one call covers both paths
    [D11,D12,D13,D22,D23,D33] = applyHessian3DAniso(I, s, spacing);
    S = sqrt(D11.^2 + D22.^2 + D33.^2 + 2*D12.^2 + 2*D13.^2 + 2*D23.^2);
    Smax = max(Smax, double(max(S(:))));
end
c = Smax / 2;
if c == 0
    c = eps;   % flat volume: avoid 0/0 in the response; R is 0 anyway
end
end
