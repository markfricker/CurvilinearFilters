function [lambda1, lambda2] = imageEigenvalues2D(I, sigma, whiteOnDark, realClass, smallEigenAbsTol)
% imageEigenvalues2D  Compute scale-normalized Hessian eigenvalues (2D)
%                      (fast masked computation only)
%
% Usage:
%   [lambda1, lambda2] = mfat.imageEigenvalues2D(I, sigma, whiteOnDark, realClass, smallEigenAbsTol)
%
% Inputs:
%   I                - HxW image (numeric), expected already or will be cast
%   sigma            - Gaussian scale (double, >0)
%   whiteOnDark      - logical; true if structures are BRIGHT on dark BG
%   realClass        - 'single' or 'double' (string) - output precision
%   smallEigenAbsTol - small-value threshold for zeroing tiny eigenvalues
%
% Outputs:
%   lambda1, lambda2 - HxW eigenvalue maps in class realClass
%                      (sorted so abs(lambda1) <= abs(lambda2))
%
% IMPLEMENTATION NOTES
% - This implementation always uses the fast mask:
%     compute inexpensive tests on Hessian (trace/determinant based),
%     find candidate pixels, compute closed-form eigenvalues only at those
%     locations, and write results back into full-size arrays.
% - This substantially reduces computation when candidate pixels are sparse
%   (typical for microscopy images with thin curvilinear structures).
% - The function is precision-aware and returns outputs cast to realClass.
%
% REFERENCES
% - S.-F. Yang and C.-H. Cheng, "Fast computation of Hessian-based
%   enhancement filters for medical images," Computer Methods and Programs
%   in Biomedicine, vol. 116, no. 3, pp. 215–225, 2014.
% - B. Obara et al., "Contrast-Independent Curvilinear Structure Detection
%   in Biomedical Images", IEEE Transactions on Image Processing, 2012.
%   (Obara-style masked eigenvalue computation inspired the optimisation.)
% - P. J. Basser, J. Mattiello, D. Le Bihan, "MR Diffusion Tensor
%   Spectroscopy and Imaging", Biophysical Journal, 1994. (fractional
%   anisotropy background)
%
% NOTE: keep the sign-convention consistent with callers:
%   if whiteOnDark == false, Hessian components are negated to match the
%   expected ridge/pivot orientation used in MFAT.
%
% -------------------------------------------------------------------------
    arguments
        I
        sigma (1,1) double {mustBePositive}
        whiteOnDark (1,1) logical = false
        realClass (1,:) char {mustBeMember(realClass,{'single','double'})} = 'single'
        smallEigenAbsTol double = 1e-4
    end

    % (1) Smooth and compute second derivatives (Hessian components)
    Iblur = gaussianSmooth2D(I, sigma);

    % central differences via gradient (keeps numeric stability)
    [Ix, Iy] = gradient(Iblur);
    [Ixx, Ixy1] = gradient(Ix);
    [Iyx, Iyy]  = gradient(Iy);
    Ixy = 0.5 .* (Ixy1 + Iyx);

    % (2) Scale-normalize Hessian
    cScale = sigma.^2;
    Hxx = cScale .* Ixx;
    Hxy = cScale .* Ixy;
    Hyy = cScale .* Iyy;

    % (3) Sign convention: flip if structures are dark on bright
    if ~whiteOnDark
        Hxx = -Hxx; Hxy = -Hxy; Hyy = -Hyy;
    end

    % Prepare output arrays (zeros by default)
    lambda1 = zeros(size(Hxx), realClass);
    lambda2 = zeros(size(Hxx), realClass);

    % ---------------------------------------------------------------------
    % fast candidate selection:
    %   B1 = - (Hxx + Hyy)   (trace-based test)
    %   B2 = Hxx.*Hyy - Hxy.^2  (determinant)
    % Select pixels where T == 1 -> plausible curvilinear candidates.
    % Compute eigenvalues only at those pixels.
    % ---------------------------------------------------------------------
    B1 = - (Hxx + Hyy);
    B2 = (Hxx .* Hyy) - (Hxy .^ 2);

    % Initial candidate map: true where potential vessel curvature exists
    T = true(size(B1));
    T(B1 < 0) = false;
    % handle degenerate case (both zero)
    T(B2 == 0 & B1 == 0) = false;

    idx = find(T);  % linear indices for candidate pixels

    if ~isempty(idx)
        % vectorize Hessian components at candidate positions
        Hxxv = Hxx(idx);
        Hxyv = Hxy(idx);
        Hyyv = Hyy(idx);

        % closed-form eigenvalues (vectorised)
        tmp = sqrt((Hxxv - Hyyv).^2 + 4 .* (Hxyv.^2));
        mu1 = 0.5 .* (Hxxv + Hyyv + tmp);
        mu2 = 0.5 .* (Hxxv + Hyyv - tmp);

        % sort by absolute value: ensure abs(lambda1) <= abs(lambda2)
        swap = abs(mu1) > abs(mu2);
        L1v = mu1;
        L2v = mu2;
        L1v(swap) = mu2(swap);
        L2v(swap) = mu1(swap);

        % cast and write back into full-size arrays
        lambda1(idx) = cast(L1v, realClass);
        lambda2(idx) = cast(L2v, realClass);
    end

    % cleanup: replace non-finite and very small values with zero
    lambda1(~isfinite(lambda1)) = cast(0, realClass);
    lambda2(~isfinite(lambda2)) = cast(0, realClass);

    lambda1(abs(lambda1) < cast(smallEigenAbsTol, realClass)) = cast(0, realClass);
    lambda2(abs(lambda2) < cast(smallEigenAbsTol, realClass)) = cast(0, realClass);

end

