function state = mfatResponseProb(cmd, response, state, opts)
% =========================================================================
% MFAT – Multiscale Fractional Anisotropy Tensor Framework
%
% Author:
%   H. Alhasson, M. Alharbi, B. Obara
%
% Refactored framework & extensions:
%   MD Fricker, Jan 2026
%
% Citation:
%   H. Alhasson, M. Alharbi, B. Obara,
%   "2D and 3D Vascular Structures Enhancement via
%    Multiscale Fractional Anisotropy Tensor",
%   ECCV Workshops (BioImage Computing), 2018.
%
% Background:
%   P. J. Basser et al.,
%   "MR Diffusion Tensor Spectroscopy and Imaging",
%   Biophysical Journal, 1994.
%
% Related work:
%   A. F. Frangi et al.,
%   "Multiscale Vessel Enhancement Filtering",
%   MICCAI, 1998.
%
% License:
%   Academic / research use. Please cite the above work.
% =========================================================================

% MFATRESPONSEPROB  Probabilistic MFAT response module
%
% OVERVIEW
%   Converts MFAT-λ responses into probabilistic evidence using Beta
%   likelihood models and log-likelihood ratio accumulation.
%
% COMMANDS
%   'init'       - Initialize probabilistic state.
%   'accumulate' - Add per-scale evidence.
%   'finalize'   - Convert LLR to posterior probability.
%
% DESIGN NOTE
%   MFAT-Prob assumes MFAT-λ has already applied strict post-processing.


switch cmd
    case 'init'
        % response size passed in
        sz = response;
        rc = opts.precision;

        state.llr = zeros(sz, rc);

        % numeric guard
        if strcmpi(rc,'single')
            state.eps0 = eps('single');
        else
            state.eps0 = eps('double');
        end

        % Beta model parameters (conservative defaults)
        % vessel: skewed toward high response
        % background: skewed toward low response
        state.vBeta = cast([2 1], rc);
        state.bBeta = cast([1 2], rc);

        state.logBV = logBetaFunc(2,1,rc);
        state.logBB = logBetaFunc(1,2,rc);

    case 'accumulate'
        r = response;

        % clamp for numerical stability
        r = min(max(r, state.eps0), 1 - state.eps0);

        % log-likelihoods
        logPv = (state.vBeta(1)-1).*log(r) + ...
                (state.vBeta(2)-1).*log(1-r) - state.logBV;

        logPb = (state.bBeta(1)-1).*log(r) + ...
                (state.bBeta(2)-1).*log(1-r) - state.logBB;

        state.llr = state.llr + (logPv - logPb);

    case 'finalize'
        % logistic transform -> posterior
        state = 1 ./ (1 + exp(-state.llr));
end
end
