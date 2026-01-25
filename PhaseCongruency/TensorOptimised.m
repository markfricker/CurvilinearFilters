function [S,Rbeta,tensor,Vx,Vy,mxc] = MDFTensor(qk,beta,c,precision)
%% MDFTensor
% Tensor-based vesselness layered on phase congruency
%
% Default precision: SINGLE
%
% Optional:
%   precision = 'single' (default) | 'double'

%% ------------------------------------------------------------
% Precision handling
%% ------------------------------------------------------------
if nargin < 4 || isempty(precision)
    precision = 'single';
end

switch lower(precision)
    case 'single'
        qk   = single(qk);
        beta = single(beta);
        c    = single(c);
        epsval = eps('single');
    case 'double'
        qk   = double(qk);
        beta = double(beta);
        c    = double(c);
        epsval = eps('double');
    otherwise
        error('Precision must be ''single'' or ''double''.');
end

%% ------------------------------------------------------------
% Dimensions and allocation
%% ------------------------------------------------------------
[m,n,theta] = size(qk);

T  = zeros(m,n,4,'like',qk);

%% ------------------------------------------------------------
% Orientation setup
%% ------------------------------------------------------------
norient = theta;
ang = ((0:norient-1)' * pi / norient);
cx  = cos(ang);
sy  = sin(ang);

%% ------------------------------------------------------------
% Tensor accumulation
%% ------------------------------------------------------------
for o = 1:theta
    q = qk(:,:,o);

    c2 = cx(o)*cx(o);
    s2 = sy(o)*sy(o);
    cs = cx(o)*sy(o);

    T(:,:,1) = T(:,:,1) + q * c2;
    T(:,:,2) = T(:,:,2) + q * cs;
    T(:,:,3) = T(:,:,3) + q * cs;
    T(:,:,4) = T(:,:,4) + q * s2;
end

%% ------------------------------------------------------------
% Analytic eigenvalues (symmetric 2x2 tensor)
%% ------------------------------------------------------------
a = T(:,:,1);
b = T(:,:,2);  % = T(:,:,3)
d = T(:,:,4);

tr    = a + d;
detT  = a.*d - b.^2;
delta = sqrt(max(tr.^2 - 4*detT, 0));

L1 = (tr - delta) / 2;
L2 = (tr + delta) / 2;

%% ------------------------------------------------------------
% Order eigenvalues so |L1| < |L2|
%% ------------------------------------------------------------
swap = abs(L1) > abs(L2);
tmp  = L1;
L1(swap) = L2(swap);
L2(swap) = tmp(swap);

%% ------------------------------------------------------------
% Eigenvector corresponding to smallest eigenvalue (vessel dir)
%% ------------------------------------------------------------
Vx = -b;
Vy = a - L1;

mag = sqrt(Vx.^2 + Vy.^2) + epsval;
Vx = Vx ./ mag;
Vy = Vy ./ mag;

%% ------------------------------------------------------------
% Vesselness metrics
%% ------------------------------------------------------------
L2(L2 == 0) = epsval;

Rbeta = L1 ./ L2;                % elongation measure
S2    = L1.^2 + L2.^2;           % tensor energy (no sqrt)

%% ------------------------------------------------------------
% Automatic scale normalization
%% ------------------------------------------------------------
Hnorm = sqrt(a.^2 + 2*b.^2 + d.^2);
mxc   = max(Hnorm(:)) / 2;

%% ------------------------------------------------------------
% Vesselness function
%% ------------------------------------------------------------
tensor = exp(-(Rbeta.^2)/(2*beta^2)) .* ...
         (1 - exp(-S2/(2*mxc^2)));

S = sqrt(S2);   % return strength for compatibility

end
