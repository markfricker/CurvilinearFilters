function options = hessian3DPresets(I, options)
%HESSIAN3DPRESETS  Safe defaults for 3D Hessian-based filters.
% 3D counterpart of hessian2DPresets.m.

if isempty(options.Sigmas)
    minDim = min([size(I,1), size(I,2), size(I,3)]);
    sigmaMax = max(2, round(minDim / 40));
    options.Sigmas = 1:1:sigmaMax;
end

if ~isfield(options,'Parameters')
    options.Parameters = struct();
end

switch lower(options.FilterType)
    case {'vesselness','plate'}
        if ~isfield(options.Parameters,'alpha'), options.Parameters.alpha = 0.5; end
        if ~isfield(options.Parameters,'beta'),  options.Parameters.beta  = 0.5; end
        % c from the data unless given (was a fixed 15 before 2026-10-03;
        % see hessian3DFrangiC). Needs the final Sigmas and Spacing.
        if ~isfield(options.Parameters,'c') || isempty(options.Parameters.c)
            spacing = [1 1 1];
            if isfield(options, 'Spacing'), spacing = options.Spacing; end
            options.Parameters.c = hessian3DFrangiC(I, options.Sigmas, spacing);
        end
end
end
