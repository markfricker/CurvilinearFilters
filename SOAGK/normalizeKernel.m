function K = normalizeKernel(K)
% Enforce zero-sum and comparable energy
K = K - mean(K(:));
K = K / (sum(abs(K(:))) + eps);
end
