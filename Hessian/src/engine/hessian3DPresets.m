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
        if ~isfield(options.Parameters,'c'),     options.Parameters.c     = 15;  end
end
end
