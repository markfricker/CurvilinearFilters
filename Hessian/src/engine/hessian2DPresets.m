function options = hessian2DPresets(I, options)
% hessian2DPresets
% Safe defaults for Hessian-based filters only

% ------------------------
% Sigma defaults
% ------------------------
if isempty(options.Sigmas)
    minDim = min(size(I,1), size(I,2));
    sigmaMax = max(2, round(minDim / 40));
    switch lower(options.FilterType)
        case 'neuriteness'
            options.Sigmas = sigmaMax;
        otherwise
            minDim = min(size(I,1), size(I,2));
            sigmaMax = max(2, round(minDim / 40));
            options.Sigmas = 1:1:sigmaMax;
    end
end

% ------------------------
% Filter-specific parameters
% ------------------------
if ~isfield(options,'Parameters')
    options.Parameters = struct();
end

switch lower(options.FilterType)

    case 'vesselness'
        if ~isfield(options.Parameters,'beta')
            options.Parameters.beta = 0.5;
        end
        if ~isfield(options.Parameters,'c')
            options.Parameters.c = 15;
        end

    case {'ridge','plate'}
        if ~isfield(options.Parameters,'alpha')
            options.Parameters.alpha = 0.5;
        end

    case {'neuriteness','blob'}
        % no parameters

end
end
