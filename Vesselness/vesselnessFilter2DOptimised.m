function [outIm, whatScale, Direction] = vesselnessFilter2DOptimised(I, options)
% This function VESSELNESSFILTER2DOPTIMISED uses the eigenvectors of the Hessian to
% compute the likeliness of an image region to vessels, according
% to the method described by Frangi:2001 (Chapter 2).
%
% [J,Scale,Direction] = vesselnessFilter2DOptimised(I, Options)
%
% inputs,
%   I : The input image (vessel image)
%   Options : Struct with input options,
%       .SigmaMin         : Minimum scale (preferred, MFAT-style)
%       .SigmaMax         : Maximum scale (preferred, MFAT-style)
%       .SigmaStep        : Scale step (preferred, MFAT-style)
%       .FrangiScaleRange : Legacy alternative to SigmaMin/SigmaMax
%       .FrangiScaleRatio : Legacy alternative to SigmaStep
%       .FrangiBetaOne    : Frangi correction constant, default 0.5
%       .FrangiBetaTwo    : Frangi correction constant, default 15
%       .WhiteOnDark      : Detect bright vessels on dark background
%       .Precision        : 'single' or 'double' (default 'single')
%       .verbose          : Show debug information, default true
%
% outputs,
%   J : The vessel enhanced image (pixel is the maximum found in all scales)
%   Scale : Matrix with the scales on which the maximum intensity 
%           of every pixel is found
%   Direction : Matrix with directions (angles) of pixels (from minor eigenvector)
%
% Reference:
%   Frangi, A. F., Niessen, W. J., Vincken, K. L., & Viergever, M. A. (1998).
%   Multiscale vessel enhancement filtering.
%   In Medical Image Computing and Computer-Assisted Intervention (MICCAI).
%
% Original Authors:
%   Marc Schrijver, 2/11/2001
%   Re-written by D. Kroon, University of Twente (May 2009)
%   Optimised implementation (streaming max, precision control), 2026

% -------------------------------------------------------------
% Default options
% -------------------------------------------------------------
defaultoptions = struct( ...
    'FrangiScaleRange', [1 10], ...
    'FrangiScaleRatio', 2, ...
    'FrangiBetaOne', 0.5, ...
    'FrangiBetaTwo', 15, ...
    'WhiteOnDark', true, ...
    'Precision', 'single', ...
    'verbose', true );

% -------------------------------------------------------------
% Process inputs
% -------------------------------------------------------------
if ~exist('options','var')
    options = defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for k = 1:numel(tags)
        if ~isfield(options,tags{k})
            options.(tags{k}) = defaultoptions.(tags{k});
        end
    end
    if numel(tags) ~= numel(fieldnames(options))
        warning('vesselnessFilter2DOptimised:unknownoption', ...
                'unknown options found');
    end
end

% % -------------------------------------------------------------
% % Align naming with MFAT-style options (backward compatible)
% % -------------------------------------------------------------
% if isfield(options,'SigmaMin') && isfield(options,'SigmaMax')
%     options.FrangiScaleRange = [options.SigmaMin options.SigmaMax];
% end
% 
% if isfield(options,'SigmaStep')
%     options.FrangiScaleRatio = options.SigmaStep;
% end
% 
% if isfield(options,'BlackWhite')
%     options.WhiteOnDark = ~options.BlackWhite;
% end

% -------------------------------------------------------------
% Precision handling
% -------------------------------------------------------------
switch lower(options.Precision)
    case 'single'
        I = single(I);
    case 'double'
        I = double(I);
    otherwise
        error('Precision must be ''single'' or ''double''');
end

% -------------------------------------------------------------
% Scale-space setup
% -------------------------------------------------------------
sigmas = options.FrangiScaleRange(1): ...
         options.FrangiScaleRatio: ...
         options.FrangiScaleRange(2);

sigmas = sort(sigmas,'ascend')
nScales = numel(sigmas);

beta = 2 * options.FrangiBetaOne^2
c    = 2 * options.FrangiBetaTwo^2
invbeta = 1 / beta;
invc    = 1 / c;

sigma2 = sigmas.^2;
rm = realmin(class(I));

% -------------------------------------------------------------
% Output buffers (streaming max, no 3D stacks)
% -------------------------------------------------------------
outIm = zeros(size(I), 'like', I);

needScale = (nargout > 1);
needDir   = (nargout > 2);

if needScale
    whatScale = zeros(size(I), 'uint16');
else
    whatScale = [];
end

if needDir
    Direction = zeros(size(I), 'like', I);
else
    Direction = [];
end

% -------------------------------------------------------------
% Vesselness filter for all sigmas
% -------------------------------------------------------------
for i = 1:nScales

    if options.verbose
        disp(['Current vesselness sigma: ' num2str(sigmas(i)) ]);
    end

    % Make 2D Hessian
    [Dxx, Dxy, Dyy] = Hessian2D(I, sigmas(i));

    % Correct for scale
    s2 = sigma2(i);
    Dxx = s2 * Dxx;
    Dxy = s2 * Dxy;
    Dyy = s2 * Dyy;

    % Calculate (abs sorted) eigenvalues and vectors
    [Lambda2, Lambda1, Ix, Iy] = eig2image(Dxx, Dxy, Dyy);

    % % Avoid division by zero
    % Lambda1s = Lambda1;
    % Lambda1s(abs(Lambda1s) < rm) = rm;

    % Compute similarity measures
    %Rb = (Lambda2 ./ Lambda1s).^2;
    Rb = (Lambda2 ./ Lambda1).^2;
    S2 = Lambda1.^2 + Lambda2.^2;

    % Vesselness response
    Ifiltered = exp(-Rb .* invbeta) .* (-expm1(-S2 .* invc));

    % Vessel polarity selection
    mask = vesselness.vesselPolarityMask(Lambda1, options.WhiteOnDark);
    Ifiltered(mask) = 0;

    % Update maximum response
    upd = Ifiltered > outIm;
    if any(upd(:))
        outIm(upd) = Ifiltered(upd);

        if needScale
            whatScale(upd) = i;
        end

        if needDir
            Direction(upd) = atan2(Ix(upd), Iy(upd));
        end
    end
end
end
