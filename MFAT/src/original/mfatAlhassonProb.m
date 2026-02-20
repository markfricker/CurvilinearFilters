function out = mfatAlhassonProb(I, sigmas, varargin)
% mfatAlhassonProb  Refined probabilistic MFAT (2D) — tuned defaults + controls
%
% Citation:
%   H. Alhasson, M. Alharbi, B. Obara,
%   "2D and 3D Vascular Structures Enhancement via
%    Multiscale Fractional Anisotropy Tensor", ECCV Workshops, 2018.
% Background: P.J. Basser et al., Biophysical Journal, 1994.
% Speedup reference: S.-F. Yang & C.-H. Cheng, Comput. Meth. Prog. Biomed., 2014.
%
% Overview:
%   This function computes a per-pixel posterior probability that a pixel
%   belongs to a curvilinear structure. For each scale:
%     - compute Hessian eigenvalues (fast Obara-style masking),
%     - form 3-factor regularised lambdas (lambda2/lambda3/lambda4),
%     - compute fractional-anisotropy-like score (faScore),
%     - map faScore -> s using same scaling as deterministic MFAT,
%     - model s with Beta distributions for vessel/background,
%     - accumulate log-likelihood ratios (LLRs) across scales,
%   then combine with prior and apply logistic to get posterior.
%
% Key improvements vs earlier draft:
%  - per-scale mapping uses same sqrt(3/2) scaling to match deterministic dynamic range
%  - added gamma and llrGain to control dynamic range robustly
%  - softened default Beta parameters to avoid extreme LLRs
%  - optional automatic polarity correction (robustLightOnDark)
%
% Usage:
%   out = mfatAlhassonProb(I, sigmas)  % sensible defaults
%   out = mfatAlhassonProb(I, sigmas, 'prior',0.1, 'gamma',0.85, 'llrGain',2.5, 'robustLightOnDark',true)
%
% Important notes on tuning:
%  - If posterior is too flat → increase llrGain (2–4) or reduce tempering (careful).
%  - If posterior is too extreme → reduce llrGain or increase gamma (toward 1).
%  - If many false positives → reduce prior or make backgroundBeta more peaked.
% -------------------------------------------------------------------------

% -----------------------
% Input parser
% -----------------------
p = inputParser;
p.FunctionName = mfilename;

addRequired(p,'I', @(x) isnumeric(x) && ismatrix(x));
addRequired(p,'sigmas', @(x) isnumeric(x) && isvector(x) && all(x>0));

addParameter(p,'tau', 0.03, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p,'tau2', 0.3, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p,'whiteOnDark', false, @(x) islogical(x) || x==0 || x==1);
addParameter(p,'prior', 0.01, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);

% softened defaults for robustness
addParameter(p,'vesselBeta', [2 1], @(x) isnumeric(x) && numel(x)==2 && all(x>0));
addParameter(p,'backgroundBeta', [1 3], @(x) isnumeric(x) && numel(x)==2 && all(x>0));

addParameter(p,'finalSmoothSigma', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p,'precision', 'single', @(x) any(validatestring(x,{'single','double'})));

% dynamic-range controls (defaults chosen conservatively)
addParameter(p,'scaleFactor', sqrt(3/2), @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'gamma', 1.0, @(x) isnumeric(x) && isscalar(x) && x>0);    % s^gamma
addParameter(p,'llrGain', 1.0, @(x) isnumeric(x) && isscalar(x) && x>0);   % multiplies tempered LLR
addParameter(p,'robustLightOnDark', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

parse(p, I, sigmas, varargin{:});

tau         = p.Results.tau;
tau2        = p.Results.tau2;
whiteOnDark = logical(p.Results.whiteOnDark);
priorProb   = p.Results.prior;
vBeta       = p.Results.vesselBeta(:)';
bBeta       = p.Results.backgroundBeta(:)';
finalSmoothSigma = p.Results.finalSmoothSigma;
precision   = validatestring(p.Results.precision, {'single','double'});

scaleFactor = p.Results.scaleFactor;
gamma       = p.Results.gamma;
llrGain     = p.Results.llrGain;
robustLightOnDark = logical(p.Results.robustLightOnDark);

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
sigmas = double(sigmas(:)');

[nRows, nCols] = size(I);
nScales = numel(sigmas);

% preallocate sum of log-likelihood ratios (LLR)
llrSum = zeros(nRows, nCols, realClass);

% precompute log Beta normalizers (cast-safe)
logBetaV = mfat.logBetaFunc(vBeta(1), vBeta(2), realClass);
logBetaB = mfat.logBetaFunc(bBeta(1), bBeta(2), realClass);

sEps = cast(1e-6, realClass);  % avoid exact 0/1 for log stability

% -----------------------
% Per-scale soft FA response and LLR accumulation
% -----------------------
for k = 1:nScales
    sigma = sigmas(k);

    % principal eigenvalue (fast masked eigenvalue computation)
    [~, lambda2] = mfat.imageEigenvalues2D(I, sigma, whiteOnDark, realClass, smallEigenAbsTol);

    % lambda3/lambda4 clipping (3-factor regularisation)
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

    % FA-like soft score (numer/denom)
    lambdaMean = abs(abs(lambda2) + abs(lambda3) + abs(lambda4)) / cast(3, realClass);

    numer = (abs(lambda2) - abs(lambdaMean)).^2 + ...
            (abs(lambda3) - abs(lambdaMean)).^2 + ...
            (abs(lambda4) - abs(lambdaMean)).^2;
    numer = sqrt(numer);

    denom = sqrt(abs(lambda2).^2 + abs(lambda3).^2 + abs(lambda4).^2);
    denom(denom == 0) = minDenom;

    faScore = numer ./ denom;
    faScore(~isfinite(faScore)) = cast(0, realClass);

    % ---------- refined mapping and contrast controls ----------
    % 1) use same scaling as deterministic MFAT to get comparable dynamic range
    s = mfat.imcomplementSafe(cast(scaleFactor, realClass) .* faScore, realClass);

    % 2) clamp to [0,1]
    s = min(max(s, cast(0, realClass)), cast(1, realClass));

    % 3) optional gamma adjustment to expand/compress midrange
    if abs(double(gamma - 1.0)) > eps
        s = s .^ cast(gamma, realClass);
    end

    % 4) clip away from exact 0/1 for log stability
    s = min(max(s, sEps), cast(1, realClass) - sEps);

    % per-pixel log-likelihoods under Beta models:
    logPv = (cast(vBeta(1) - 1, realClass)) .* log(s) + (cast(vBeta(2) - 1, realClass)) .* log(1 - s) - logBetaV;
    logPb = (cast(bBeta(1) - 1, realClass)) .* log(s) + (cast(bBeta(2) - 1, realClass)) .* log(1 - s) - logBetaB;

    llr = logPv - logPb;

    % accumulate evidence across scales
    llrSum = llrSum + llr;
end

% -----------------------
% Temper the summed LLR to avoid overconfident aggregation
% -----------------------
if nScales > 1
    llrSum = llrSum ./ cast(sqrt(double(nScales)), realClass);
end

% -----------------------
% Apply llrGain to adjust dynamic range in log-odds
% -----------------------
if abs(double(llrGain - 1.0)) > eps
    llrSum = llrSum .* cast(llrGain, realClass);
end

% -----------------------
% Combine with prior to form posterior probability
% -----------------------
priorProb = min(max(priorProb, tinyEps), 1 - tinyEps);
logPriorOdds = log(double(priorProb) ./ (1 - double(priorProb)));
logPriorOdds = cast(logPriorOdds, realClass);

postLogOdds = logPriorOdds + llrSum;

% stable logistic transform
post = mfat.logisticStable(postLogOdds, realClass);

% optional final smoothing
if finalSmoothSigma > 0
    try
        post = imgaussfilt(post, finalSmoothSigma);
    catch
        h = fspecial('gaussian', max(3,2*ceil(3*finalSmoothSigma)+1), finalSmoothSigma);
        post = imfilter(post, h, 'same', 'replicate');
    end
    post = min(max(post, cast(0, realClass)), cast(1, realClass));
end

% -----------------------
% Optional robust polarity check and flip (median-based)
% -----------------------
if robustLightOnDark
    % choose a high percentile (top 0.5%) as candidate vessel pixels
    pctileVal = 99.5;
    pThr = prctile(post(:), pctileVal);
    vesselCand = post >= pThr;

    minFrac = 5e-4; % 0.05% of image
    if nnz(vesselCand) >= max(3, round(minFrac * numel(post)))
        medV = median(I(vesselCand));          % median intensity on candidate pixels
        medBg = median(I(~vesselCand));       % median intensity on remainder

        % if disagreement with whiteOnDark -> flip
        if (medV < medBg && whiteOnDark) || (medV > medBg && ~whiteOnDark)
            post = 1 - post;
        end
    end
end

out = post;
end
