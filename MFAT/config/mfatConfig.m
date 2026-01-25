function cfg = mfatConfig()
% MFATCONFIG  Default configuration for MFAT framework
%
% OVERVIEW
%   Returns a struct containing default parameters for MFAT-λ, MFAT-Prob,
%   entropy weighting, and fractional shaping.
%
% DESIGN
%   Centralizes all parameters to ensure reproducibility and clarity.

cfg.core.tau         = 0.03;
cfg.core.tau2        = 0.3;
cfg.core.D           = 0.27;
cfg.core.whiteOnDark = true;
cfg.core.precision   = 'single';

cfg.prob.enabled     = true;
cfg.prob.vBeta       = [2 1];
cfg.prob.bBeta       = [1 2];
cfg.prob.llrGain     = 1.0;
cfg.prob.prior       = 0.5;

cfg.entropy.enabled  = false;
cfg.entropy.beta     = 0.5;
cfg.entropy.combine  = 'mean';

cfg.fractional.enabled = false;
cfg.fractional.alpha   = 0.7;
end
