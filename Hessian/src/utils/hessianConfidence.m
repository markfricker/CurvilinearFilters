function C = hessianConfidence(resp, l1, l2)
% hessianConfidence
% Confidence for Hessian-based responses

resp = single(resp);
epsv = eps('single');

% 1) response strength
Cr = resp / (max(resp(:)) + epsv);

% 2) anisotropy confidence
Ca = abs(l2) ./ (abs(l1) + abs(l2) + epsv);

% 3) combine conservatively
C = Cr .* Ca;

% clamp
C(C < 0) = 0;
C(C > 1) = 1;
end
