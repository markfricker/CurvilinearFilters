function out = mfatAlhassonLambda(I, sigmas, varargin)
% mfatAlhassonLambda  Deterministic 3-factor MFAT (2D) - lambda version
%
% -------------------------------------------------------------------------
% @cite
% -------------------------------------------------------------------------
% If you use this code in academic work, please cite:
% [1] H. Alhasson, M. Alharbi, B. Obara,
%     "2D and 3D Vascular Structures Enhancement via
%      Multiscale Fractional Anisotropy Tensor",
%     ECCV Workshops (BioImage Computing), 2018.
%
% Origin of fractional anisotropy:
% [2] P. J. Basser, J. Mattiello, D. Le Bihan,
%     "MR Diffusion Tensor Spectroscopy and Imaging",
%     Biophysical Journal, 1994.
%
% Foundational vessel enhancement:
% [3] A. F. Frangi, W. J. Niessen, K. L. Vincken, M. A. Viergever,
%     "Multiscale Vessel Enhancement Filtering", MICCAI, 1998.
%
% Additional related work:
% [4] B. Obara, M. Fricker, D. Gavaghan, V. Grau,
%     "Contrast-Independent Curvilinear Structure Detection in Biomedical Images",
%     IEEE Trans. Image Processing, 2012.
% -------------------------------------------------------------------------
%
% Overview
% -------------------------------------------------------------------------
% mfatAlhassonLambda implements the deterministic Multiscale Fractional
% Anisotropy Tensor (MFAT) filter using the original three-factor eigenvalue
% regularisation. For each Gaussian scale sigma the code:
%   1) computes the Hessian and scale-normalized second derivatives;
%   2) extracts principal eigenvalue lambda2 and forms regularised lambda3
%      and lambda4 using tau and tau2 thresholds (3-factor);
%   3) computes a fractional-anisotropy-like response from lambda2/lambda3/lambda4;
%   4) inverts & applies strict Obara/Alhasson post-processing rules;
%   5) aggregates across scales with a deterministic D * tanh(...) update,
%      combined with per-scale clamping to ensure stability and reproducibility.
%
% This deterministic formulation reproduces the predictable behaviour useful
% for dense networks (e.g., ER) where binary-like responses and reproducibility
% are desired.
% -------------------------------------------------------------------------
%
% Parameters, meaning & guidance
% -------------------------------------------------------------------------
%  - sigmas (vector): Gaussian scales to probe. Pick a range covering about
%       0.5× to 2× the expected object half-width.
%       Examples:
%         thin tubules (ER):     0.5 : 0.5 : 3
%         medium curvilinear:    1   : 1   : 6
%         blobs/vesicles:        2   : 1   : 8
%
%  - tau (scalar, default 0.03): lower clipping factor for lambda3.
%       Smaller tau => more permissive to thin structures (try 0.005–0.03).
%  - tau2 (scalar, default 0.3): larger clipping factor for lambda4; controls
%       tolerance to blob-like curvature (try 0.1–0.6).
%  - D (scalar, default 0.27): multiscale aggregation step. Larger D ->
%       stronger accumulation and crisper/binary outputs. Typical: 0.15–0.4.
%  - whiteOnDark (logical): true if structures are BRIGHT on dark BG, false
%       if DARK on bright BG (default false). This flips Hessian sign.
%  - precision ('single'|'double'): default 'single'. 'single' recommended for
%       speed/memory; use 'double' for debugging if needed.
%
% Tuning recipes:
%  - emphasize thin tubes: sigmas = 0.5:0.5:3, tau ≈ 0.005–0.03, tau2 ≈ 0.1–0.3, D ≈ 0.20–0.30
%  - shift to oblong mitochondria: sigmas = 1:0.5:4, tau ≈ 0.03–0.08, tau2 ≈ 0.2–0.5, D ≈ 0.25–0.35
%  - emphasize round vesicles: include larger sigmas up to 6–8, tau ≈ 0.05–0.2, tau2 ≈ 0.4–0.8, increase D for crispness
%
% Troubleshooting:
%  - output inverted? Check whiteOnDark setting.
%  - output too binary? Reduce D or apply mild final gaussian smoothing externally.
%  - missing thin edges? Decrease tau and include smaller sigmas.
% -------------------------------------------------------------------------
%
% USAGE:
%  out = mfatAlhassonLambda(I, sigmas)
%  out = mfatAlhassonLambda(I, sigmas, 'tau',0.03,'tau2',0.3,'D',0.27,'whiteOnDark',false,'precision','single')
%
% Inputs:
%  I      - 2D grayscale image (numeric)
%  sigmas - vector of Gaussian scales (double)
%
% Output:
%  out - vesselness map in chosen precision, normalized to [0,1]
% -------------------------------------------------------------------------

% -----------------------
% Input parsing
% -----------------------
p = inputParser;
p.FunctionName = mfilename;
addRequired(p,'I', @(x) isnumeric(x) && ismatrix(x));
addRequired(p,'sigmas', @(x) isnumeric(x) && isvector(x) && all(x>0));
addParameter(p,'tau', 0.03, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p,'tau2', 0.3, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p,'D', 0.27, @(x) isnumeric(x) && isscalar(x));
addParameter(p,'whiteOnDark', true, @(x) islogical(x) || x==0 || x==1);
addParameter(p,'precision', 'single', @(x) any(validatestring(x,{'single','double'})));
parse(p, I, sigmas, varargin{:});

tau         = p.Results.tau;
tau2        = p.Results.tau2;
D           = p.Results.D;
whiteOnDark = logical(p.Results.whiteOnDark);
precision   = validatestring(p.Results.precision, {'single','double'});

% -----------------------
% Precision & numeric guards
% -----------------------
if strcmpi(precision,'single')
    realClass = 'single';
    tinyEps = eps('single');
    smallEigenAbsTol = single(1e-4);
else
    realClass = 'double';
    tinyEps = eps('double');
    smallEigenAbsTol = 1e-8;
end
minDenom = cast(tinyEps, realClass);

% -----------------------
% Preprocess image (cast & normalize)
% -----------------------
I = cast(I, realClass);
maxI = max(I(:));
if maxI > 0
    I = I ./ cast(maxI, realClass);
end
sigmas = double(sigmas(:)');  % allow fractional sigmas

% -----------------------
% Core MFAT (lambda-driven deterministic update)
% -----------------------
vesselness = zeros(size(I), realClass);

for idxScale = 1:numel(sigmas)
    sigma = sigmas(idxScale);

    % (1) principal Hessian eigenvalue (lambda2)
    [~, lambda2] = mfat.imageEigenvalues2D(I, sigma, whiteOnDark, realClass, smallEigenAbsTol);

    % (2) regularised eigenvalues (lambda3, lambda4)
    lambda3 = lambda2;
    min3 = min(lambda3(:));
    if min3 < 0
        mask3 = (lambda3 < 0) & (lambda3 >= cast(tau * min3, realClass));
        lambda3(mask3) = cast(tau * min3, realClass);
    end

    lambda4 = lambda2;
    min4 = min(lambda4(:));
    if min4 < 0
        mask4 = (lambda4 < 0) & (lambda4 >= cast(tau2 * min4, realClass));
        lambda4(mask4) = cast(tau2 * min4, realClass);
    end

    % (3) fractional-anisotropy-like response
    lambdaMean = abs(abs(lambda2) + abs(lambda3) + abs(lambda4)) / cast(3, realClass);
    numer = (abs(lambda2) - abs(lambdaMean)).^2 + ...
            (abs(lambda3) - abs(lambdaMean)).^2 + ...
            (abs(lambda4) - abs(lambdaMean)).^2;
    numer = sqrt(numer);

    denom = sqrt(abs(lambda2).^2 + abs(lambda3).^2 + abs(lambda4).^2);
    denom(denom == 0) = minDenom;

    response = numer ./ denom;
    response(~isfinite(response)) = cast(0, realClass);

    % (4) inversion / scaling
    response = mfat.imcomplementSafe(cast(sqrt(cast(3, realClass)./cast(2, realClass)), realClass) .* response, realClass);
    response(~isfinite(response)) = cast(0, realClass);

    % (5) strict post-processing (Obara/Alhasson rules)
    x = lambda3 - lambda2;
    minX = min(x(:));
    maxX = max(x(:));
    response(x == minX) = cast(1, realClass);
    response(x < maxX)    = cast(0, realClass);
    response(lambda2 > x) = cast(0, realClass);
    response(lambda3 > x) = cast(0, realClass);
    response(lambda2 >= 0) = cast(0, realClass);
    response(lambda3 >= 0) = cast(0, realClass);
    response(~isfinite(response)) = cast(0, realClass);

    % (6) multiscale update (D * tanh) and per-scale clamp
    if idxScale == 1
        vesselness = response;
    else
        vesselness = vesselness + cast(D, realClass) .* tanh(response - cast(D, realClass));
        vesselness = max(vesselness, response);
    end

    vesselness = min(max(vesselness, cast(0, realClass)), cast(1, realClass));
end

% -----------------------
% Final normalization & cleanup
% -----------------------
out = vesselness;
maxOut = max(out(:));
if maxOut > 0
    out = out ./ cast(maxOut, realClass);
end
out(out < cast(1e-2, realClass)) = cast(0, realClass);

end
