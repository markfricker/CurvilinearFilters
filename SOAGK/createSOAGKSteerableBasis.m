function [Gxx, Gyy, Gxy] = createSOAGKSteerableBasis(sigma, params, precision)

% --- Image precision ---
castfun = str2func(precision);

n = params.size;
h = (n-1)/2;
[x, y] = meshgrid(-h:h, -h:h);

x = castfun(x);
y = castfun(y);

rho  = castfun(params.rho);
irho = 1/rho;

sigma2 = castfun(sigma^2);
sigma4 = sigma2^2;

% --- Anisotropic coordinates ---
x2 = rho  * x;
y2 = irho * y;

r2 = x2.^2 + y2.^2;
G  = exp(-r2/(2*sigma2));

Gxx = (x2.^2/sigma4 - 1/sigma2) .* G;
Gyy = (y2.^2/sigma4 - 1/sigma2) .* G;
Gxy = (x2.*y2/sigma4) .* G;

% --- Normalize ---
Gxx = normalizeKernel(Gxx);
Gyy = normalizeKernel(Gyy);
Gxy = normalizeKernel(Gxy);

% 🔴 CRITICAL FIX: cast kernels to double
Gxx = double(Gxx);
Gyy = double(Gyy);
Gxy = double(Gxy);
end

