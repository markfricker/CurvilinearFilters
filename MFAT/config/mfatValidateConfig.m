function mfatValidateConfig(cfg)

% ---- core ----
assert(cfg.core.tau >= 0, 'tau must be >= 0');
assert(cfg.core.tau2 >= 0, 'tau2 must be >= 0');
assert(isfinite(cfg.core.D), 'D must be finite');
assert(ismember(cfg.core.precision, {'single','double'}), ...
    'precision must be single or double');

% ---- prob ----
assert(numel(cfg.prob.vBeta)==2 && all(cfg.prob.vBeta>0), ...
    'vBeta must be [a b] > 0');
assert(numel(cfg.prob.bBeta)==2 && all(cfg.prob.bBeta>0), ...
    'bBeta must be [a b] > 0');
assert(cfg.prob.llrGain > 0, 'llrGain must be > 0');
assert(cfg.prob.prior > 0 && cfg.prob.prior < 1, ...
    'prior must be in (0,1)');

% ---- entropy ----
assert(cfg.entropy.beta >= 0, 'entropy beta must be >= 0');
assert(ismember(cfg.entropy.combine, {'mean','min','max'}), ...
    'entropy.combine must be mean|min|max');

% ---- fractional ----
assert(cfg.fractional.alpha > 0, 'fractional alpha must be > 0');

end
