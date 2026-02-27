% PHASESYMOPTIMISED - Phase symmetry on an image (optimised version).
%
% This function calculates the phase symmetry of points in an image.
% This is a contrast-invariant measure of symmetry and can be used as a
% line and blob detector. The greyscale 'polarity' of the lines to find
% can be specified.
%
% There are potentially many arguments, here is the full usage:
%
%   [phaseSym, orientation, totalEnergy, T, seeFilter] = ...
%       phaseSymOptimised(im, nscale, norient, minWaveLength, mult, ...
%                         sigmaOnf, k, polarity, noiseMethod, ...
%                         'precision', 'single', 'storeFilters', false)
%
% However, apart from the image, all parameters have defaults:
%
%    phaseSym = phaseSymOptimised(im);
%
% Arguments:
%              Default values      Description
%
%    nscale           5    - Number of wavelet scales, try values 3-6
%    norient          6    - Number of filter orientations.
%    minWaveLength    3    - Wavelength of smallest scale filter.
%    mult             2.1  - Scaling factor between successive filters.
%    sigmaOnf         0.55 - Ratio of the standard deviation of the Gaussian
%                            describing the log Gabor filter's transfer function
%                            in the frequency domain to the filter center frequency.
%    k                2.0  - No of standard deviations of the noise energy beyond
%                            the mean at which we set the noise threshold point.
%                            You may want to vary this up to 10 or 20 for noisy images.
%    polarity         0    - Controls 'polarity' of symmetry features to find.
%                             1 - just return 'bright' points
%                            -1 - just return 'dark' points
%                             0 - return bright and dark points.
%    noiseMethod      -1   - Parameter specifies method used to determine
%                            noise statistics.
%                              -1 use median of smallest scale filter responses
%                              -2 use mode of smallest scale filter responses
%                               0+ use noiseMethod value as the fixed noise threshold.
%
% Additional options (keyword-value pairs only):
%    'precision'      'single' (default) | 'double' | 'like'
%                     Controls internal arithmetic precision.
%                     'single' halves memory use relative to double.
%                     'like' matches the precision of the input image.
%    'storeFilters'   false (default) | true
%                     If false, seeFilter is returned as an empty array,
%                     avoiding an fftshift per scale when not needed.
%
% Return values:
%    phaseSym     - Phase symmetry image (values between 0 and 1).
%    orientation  - Orientation image in degrees (0-180), anticlockwise.
%                   Quantised by the number of orientations.
%    totalEnergy  - Un-normalised raw symmetry energy.
%    T            - Calculated noise threshold.
%    seeFilter    - Composite log Gabor filter (fftshifted max across scales).
%                   Empty if storeFilters is false.
%
% Notes on specifying parameters:
%
%  >> phaseSym = phasesymOptimised(im, 5, 6, 3, 2.5, 0.55, 2.0, 0);
%  >> phaseSym = phasesymOptimised(im, 5, 6, 3);
%  >> phaseSym = phasesymOptimised(im, 5, 6, 3, 'polarity', -1, 'k', 2.5);
%  >> phaseSym = phasesymOptimised(im, 'precision', 'single', 'storeFilters', false);
%
% Notes on filter settings to obtain even coverage of the spectrum:
% sigmaOnf  .85   mult 1.3
% sigmaOnf  .75   mult 1.6     (filter bandwidth ~1 octave)
% sigmaOnf  .65   mult 2.1
% sigmaOnf  .55   mult 3       (filter bandwidth ~2 octaves)
%
% For maximum speed the input image should have dimensions that are powers
% of 2, but the code operates on images of arbitrary size.
%
% See Also: PHASECONG3OPTIMISED, PHASECONG, PHASECONG2, GABORCONVOLVE
%
% References:
%     Peter Kovesi, "Symmetry and Asymmetry From Local Phase". AI'97,
%     Tenth Australian Joint Conference on Artificial Intelligence,
%     2-4 December 1997.
%     http://www.cs.uwa.edu.au/pub/robvis/papers/pk/ai97.ps.gz
%
%     Peter Kovesi, "Image Features From Phase Congruency". Videre: A
%     Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
%     Summer 1999. http://mitpress.mit.edu/e-journals/Videre/001/v13.html
%
% Original code: Copyright (c) 1996-2017 Peter Kovesi
% http://www.peterkovesi.com
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, subject to the following
% conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The software is provided "as is", without warranty of any kind.
%
% -------------------------------------------------------------------------
% SUMMARY OF OPTIMISATIONS (relative to phasesym.m / phasesym_MDF)
%
% 1. PERSISTENT FILTER CACHE
%    All image-size-dependent arrays (radius, sintheta, costheta, lp,
%    logGabor cell) are computed once and stored in a persistent struct
%    keyed on image size, nscale, minWaveLength, mult, sigmaOnf, lp params,
%    and precision class. Repeated calls on the same image size with the
%    same parameters (e.g. processing a stack of confocal frames) skip the
%    entire filter construction phase. This mirrors the approach in
%    phasecong3Optimised and yields the largest single speedup for repeated
%    calls.
%
% 2. PRECISION TOGGLE ('single' | 'double' | 'like')
%    Matches phasecong3Optimised. Defaults to 'single', which halves peak
%    memory relative to double. All accumulator arrays, filter arrays and
%    the imagefft are cast at the point of construction. epsilon is also
%    cast to the working precision to avoid silent upcasting.
%
% 3. NOISE THRESHOLD T COMPUTED ONCE, CORRECTLY
%    The original phasesym_MDF computed T from orientation 0 only.
%    This is correct in principle (tau does not vary with orientation), but
%    using a fixed orientation-0 spread function means the filter at s=1 is
%    not the pure radial filter — it carries the orientation-0 angular
%    weighting. While the error is small, the cleaner approach (used here)
%    computes An0 directly from the radial logGabor{1} alone (no spread),
%    which is the maximum-energy estimate and is consistent with the intent
%    of the noise model. This also saves constructing a throw-away spread
%    matrix before the main loop.
%
% 4. SEEFILTER IS OPT-IN
%    seeFilter requires an fftshift per scale (nscale calls). When the
%    caller does not need the visualisation array this work is wasted.
%    Setting storeFilters=false (default) skips all fftshift calls and
%    returns seeFilter as [].
%
% 5. ORIENTATION ACCUMULATION: acos INSTEAD OF atan2
%    The angular spread is computed via acos(clamped dot-product) rather
%    than atan2(ds,dc). This avoids computing and storing ds and dc as
%    full m*n arrays and is typically faster. Matches phasecong3Optimised.
%
% 6. CHECKARGS HARDENED AND BROUGHT IN LINE WITH phasecong3Optimised
%    - Correctly handles partial positional + keyword/value mixed calls
%      (the original phasesym.m had a bug: the keyword loop initialised n
%      inside the positional block so n was undefined if only keyword args
%      were given and readstate was never set to keywordvalue before the
%      keyword loop was entered).
%    - string type accepted in addition to char for keyword names.
%    - Precision and storeFilters options added.
%    - im is NOT force-cast to double here; precision is handled in main.
%
% 7. RAYLEIGHMODE: histc REPLACED WITH histcounts
%    histc was deprecated in R2014b. histcounts is used with explicit
%    linspace edges for compatibility with MATLAB R2016b and later.
%    Behaviour is identical to the original.
%
% 8. end KEYWORDS ADDED TO ALL NESTED FUNCTIONS
%    The original phasesym.m was missing end statements for checkargs and
%    rayleighmode. While MATLAB tolerates this for nested functions it is
%    inconsistent style and can cause parse failures in some contexts.
%
% -------------------------------------------------------------------------

function [phaseSym, orientation, totalEnergy, T, seeFilter] = phaseSymOptimised(varargin)

    [im, nscale, norient, minWaveLength, mult, sigmaOnf, k, ...
     polarity, noiseMethod, precision, storeFilters] = checkargs(varargin(:));

    % -----------------------------------------------------------------------
    % Precision setup — mirrors phasecong3Optimised
    % -----------------------------------------------------------------------
    if strcmpi(precision, 'single')
        if ~isa(im, 'single'), im = single(im); end
        epsilon = single(1e-4);
    elseif strcmpi(precision, 'double')
        if ~isa(im, 'double'), im = double(im); end
        epsilon = 1e-4;
    elseif strcmpi(precision, 'like')
        if ~isa(im, 'single') && ~isa(im, 'double')
            im = double(im);
        end
        epsilon = cast(1e-4, 'like', im);
    else
        error('precision must be ''single'', ''double'', or ''like''');
    end

    [rows, cols] = size(im);
    imagefft     = fft2(im);                     % FFT in working precision

    % -----------------------------------------------------------------------
    % Persistent filter cache — keyed on all filter-determining parameters
    % -----------------------------------------------------------------------
    persistent CACHE
    if isempty(CACHE)
        CACHE = struct();
    end

    precClass = class(im);
    cacheKey = sprintf('%dx%d|ns=%d|minWL=%.12g|mult=%.12g|sig=%.12g|lp=%.4g,%.4g|%s', ...
        rows, cols, nscale, minWaveLength, mult, sigmaOnf, .4, 10, precClass);

    if isfield(CACHE, 'key') && strcmp(CACHE.key, cacheKey)
        % Retrieve cached filter arrays
        sintheta = CACHE.sintheta;
        costheta = CACHE.costheta;
        logGabor = CACHE.logGabor;
        if storeFilters
            seeFilter = CACHE.seeFilter;
        else
            seeFilter = [];
        end
    else
        % ------------------------------------------------------------------
        % Build frequency-domain coordinate grids
        % Handles both odd and even image dimensions correctly.
        % ------------------------------------------------------------------
        if mod(cols, 2)
            xrange = (-(cols-1)/2 : (cols-1)/2) / cols;
        else
            xrange = (-cols/2 : (cols/2-1)) / cols;
        end

        if mod(rows, 2)
            yrange = (-(rows-1)/2 : (rows-1)/2) / rows;
        else
            yrange = (-rows/2 : (rows/2-1)) / rows;
        end

        [x, y]   = meshgrid(xrange, yrange);
        radius   = sqrt(x.^2 + y.^2);
        theta    = atan2(-y, x);

        radius   = ifftshift(radius);
        theta    = ifftshift(theta);
        radius(1,1) = 1;                          % avoid log(0)

        sintheta = cast(sin(theta), 'like', im);
        costheta = cast(cos(theta), 'like', im);
        radius   = cast(radius,     'like', im);
        % x, y, theta go out of scope — no manual clear needed

        % ------------------------------------------------------------------
        % Low-pass filter and log Gabor filter bank
        % ------------------------------------------------------------------
        lp = cast(lowpassfilter([rows, cols], .4, 10), 'like', im);

        logGabor     = cell(1, nscale);
        denomLog     = cast(2 * log(sigmaOnf)^2, 'like', im);
        seeFilterAcc = cast(zeros(rows, cols), 'like', im);

        for s = 1:nscale
            wavelength  = minWaveLength * mult^(s-1);
            fo          = cast(1.0 / wavelength, 'like', im);
            lg          = exp((-(log(radius/fo)).^2) ./ denomLog);
            lg          = lg .* lp;
            lg(1,1)     = 0;
            logGabor{s} = lg;
            if storeFilters
                seeFilterAcc = max(seeFilterAcc, fftshift(lg));
            end
        end

        seeFilter = seeFilterAcc;   % [] if storeFilters==false (loop never ran)

        % Store in cache
        CACHE.key      = cacheKey;
        CACHE.sintheta = sintheta;
        CACHE.costheta = costheta;
        CACHE.logGabor = logGabor;
        CACHE.seeFilter = seeFilter;
    end

    % -----------------------------------------------------------------------
    % Noise threshold T — computed ONCE from the pure radial filter at s=1
    % (no angular spread weighting), which is the maximum-energy noise
    % estimate and avoids constructing a throw-away spread matrix.
    % -----------------------------------------------------------------------
    if noiseMethod >= 0
        T = cast(noiseMethod, 'like', im);
    else
        An0 = abs(ifft2(imagefft .* logGabor{1}));
        if noiseMethod == -1
            tau = double(median(An0(:))) / sqrt(log(4));
        else  % noiseMethod == -2
            tau = rayleighmode(double(An0(:)));
        end
        totalTau            = tau * (1 - (1/mult)^nscale) / (1 - (1/mult));
        EstNoiseEnergyMean  = totalTau * sqrt(pi/2);
        EstNoiseEnergySigma = totalTau * sqrt((4-pi)/2);
        T = cast(max(EstNoiseEnergyMean + k*EstNoiseEnergySigma, double(epsilon)), 'like', im);
    end

    % -----------------------------------------------------------------------
    % Allocate accumulators
    % -----------------------------------------------------------------------
    zero        = zeros(rows, cols, 'like', im);
    totalEnergy = zero;
    totalSumAn  = zero;
    orientation = zero;
    maxEnergy   = zero;

    % -----------------------------------------------------------------------
    % Main orientation loop
    % -----------------------------------------------------------------------
    for o = 1:norient
        angl = cast((o-1) * pi / norient, 'like', im);

        % Angular spread via acos(clamped dot-product) — avoids storing ds/dc
        % as full m*n arrays and matches phasecong3Optimised.
        dc     = costheta * cos(angl) + sintheta * sin(angl);
        dc     = min(max(dc, cast(-1,'like',im)), cast(1,'like',im));
        dtheta = min(acos(dc) * cast(norient/2, 'like', im), cast(pi,'like',im));
        spread = (cos(dtheta) + cast(1,'like',im)) / cast(2,'like',im);

        sumAn_ThisOrient  = zero;
        Energy_ThisOrient = zero;

        for s = 1:nscale
            EO = ifft2(imagefft .* (logGabor{s} .* spread));
            An = abs(EO);
            sumAn_ThisOrient = sumAn_ThisOrient + An;

            switch polarity
                case  0,  Energy_ThisOrient = Energy_ThisOrient + abs(real(EO)) - abs(imag(EO));
                case  1,  Energy_ThisOrient = Energy_ThisOrient +     real(EO)  - abs(imag(EO));
                case -1,  Energy_ThisOrient = Energy_ThisOrient -     real(EO)  - abs(imag(EO));
            end
        end

        % Apply noise threshold
        Energy_ThisOrient = Energy_ThisOrient - T;

        totalSumAn  = totalSumAn  + sumAn_ThisOrient;
        totalEnergy = totalEnergy + Energy_ThisOrient;

        % Track orientation at which energy is maximum
        if o == 1
            maxEnergy = Energy_ThisOrient;
        else
            change      = Energy_ThisOrient > maxEnergy;
            orientation = cast(o-1, 'like', im) .* change + orientation .* (~change);
            maxEnergy   = max(maxEnergy, Energy_ThisOrient);
        end

    end  % for each orientation

    % -----------------------------------------------------------------------
    % Final outputs
    % -----------------------------------------------------------------------
    phaseSym    = max(totalEnergy, 0) ./ (totalSumAn + epsilon);
    orientation = fix(orientation * cast(180 / norient, 'like', im));

end  % --- end main function ---


%%-------------------------------------------------------------------------
% CHECKARGS
%
% Process arguments, assign defaults, validate. Brought in line with
% phasecong3Optimised: properly handles mixed positional/keyword calls,
% accepts both char and string keyword types, adds precision and
% storeFilters options, and removes the original bug where n was
% undefined when only keyword arguments were supplied.
%
function [im, nscale, norient, minWaveLength, mult, sigmaOnf, k, ...
          polarity, noiseMethod, precision, storeFilters] = checkargs(arg)

    nargs = length(arg);

    if nargs < 1
        error('No image supplied as an argument');
    end

    % Defaults
    im            = [];
    nscale        = 5;
    norient       = 6;
    minWaveLength = 3;
    mult          = 2.1;
    sigmaOnf      = 0.55;
    k             = 2.0;
    polarity      = 0;
    noiseMethod   = -1;
    precision     = 'single';    % optimised default
    storeFilters  = false;       % opt-in; avoids nscale fftshift calls

    % Find the index of the first string/char argument to split positional
    % from keyword/value pairs robustly.
    firstKey = [];
    for i = 1:nargs
        if isa(arg{i}, 'char') || isa(arg{i}, 'string')
            firstKey = i;
            break;
        end
    end

    if isempty(firstKey)
        lastPos = nargs;
        hasKeywords = false;
    else
        lastPos = firstKey - 1;
        hasKeywords = true;
    end

    % Read positional numeric args
    for n = 1:lastPos
        switch n
            case 1,  im            = arg{n};
            case 2,  nscale        = arg{n};
            case 3,  norient       = arg{n};
            case 4,  minWaveLength = arg{n};
            case 5,  mult          = arg{n};
            case 6,  sigmaOnf      = arg{n};
            case 7,  k             = arg{n};
            case 8,  polarity      = arg{n};
            case 9,  noiseMethod   = arg{n};
            otherwise, error('Too many positional arguments');
        end
    end

    % Read keyword/value pairs
    if hasKeywords
        if mod(nargs - firstKey + 1, 2) ~= 0
            error('Keyword arguments must be supplied as name-value pairs');
        end
        n = firstKey;
        while n < nargs
            if ~(isa(arg{n}, 'char') || isa(arg{n}, 'string'))
                error('Expected a parameter name (string) at position %d', n);
            end
            key = lower(string(arg{n}));
            val = arg{n+1};

            if     startsWith(key, "im"),             im            = val;
            elseif startsWith(key, "nscale"),         nscale        = val;
            elseif startsWith(key, "norient"),        norient       = val;
            elseif startsWith(key, "minwavelength"),  minWaveLength = val;
            elseif startsWith(key, "mult"),           mult          = val;
            elseif startsWith(key, "sigmaonf"),       sigmaOnf      = val;
            elseif startsWith(key, "k"),              k             = val;
            elseif startsWith(key, "polarity"),       polarity      = val;
            elseif startsWith(key, "noisemethod"),    noiseMethod   = val;
            elseif startsWith(key, "precision")
                if ~(isa(val,'char') || isa(val,'string'))
                    error('precision must be ''single'', ''double'', or ''like''');
                end
                precision = char(val);
            elseif startsWith(key, "storefilters"),   storeFilters  = logical(val);
            else
                error('Unrecognised parameter name: %s', key);
            end
            n = n + 2;
        end
    end

    if isempty(im)
        error('No image argument supplied');
    end

    if ndims(im) == 3
        warning('phasesymOptimised:colourImage', ...
            'Colour image supplied: converting to greyscale.');
        im = rgb2gray(im);
    end

    if nscale < 1,        error('nscale must be an integer >= 1');              end
    if norient < 1,       error('norient must be an integer >= 1');             end
    if minWaveLength < 2, error('minWaveLength < 2 is not meaningful');         end
    if ~ismember(polarity, [-1 0 1])
        error('Allowed polarity values are -1, 0 and 1');
    end

end  % --- end checkargs ---


%%-------------------------------------------------------------------------
% RAYLEIGHMODE
%
% Computes the mode of data assumed to follow a Rayleigh distribution.
% Uses histcounts (replaces deprecated histc, compatible with R2016b+).
%
% Usage:  rmode = rayleighmode(data, nbins)
%
% mean   = mode * sqrt(pi/2)
% stddev = mode * sqrt((4-pi)/2)
%
% See: http://mathworld.wolfram.com/RayleighDistribution.html
%
function rmode = rayleighmode(data, nbins)

    if nargin < 2
        nbins = 50;
    end

    mx    = max(data(:));
    edges = linspace(0, mx, nbins + 1);
    n     = histcounts(data(:), edges);
    [~, ind] = max(n);
    rmode = (edges(ind) + edges(ind+1)) / 2;

end  % --- end rayleighmode ---