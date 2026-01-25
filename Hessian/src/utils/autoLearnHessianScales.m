function sigmas = autoLearnHessianScales(I, sigmaMin, sigmaMax, sigmaStep)
% autoLearnHessianScales
% Selects informative scales based on energy change

I = single(I);
sigmasAll = sigmaMin:sigmaStep:sigmaMax;

energy = zeros(size(sigmasAll),'single');

for k = 1:numel(sigmasAll)
    s = sigmasAll(k);

    [Dxx,~,Dyy] = Hessian2D(I, s);
    curv = abs(s^2 * Dxx) + abs(s^2 * Dyy);

    energy(k) = mean(curv(:));
end

% Normalise
energy = energy / max(energy);

% Compute relative change
dE = diff(energy) ./ energy(1:end-1);

% Always keep smallest scale
keep = false(size(sigmasAll));
keep(1) = true;

% Keep scales with significant incremental contribution
keep(2:end) = dE > 0.05;   % 5% relative increase

% Safety: ensure at least one scale
if ~any(keep)
    keep(1) = true;
end

sigmas = sigmasAll(keep);
end



