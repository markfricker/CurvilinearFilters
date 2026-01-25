[xs, ys] = meshgrid(-nHalfSize:nHalfSize, -nHalfSize:nHalfSize);

ct = cos(params.theta);
st = sin(params.theta);
rho = params.rho;
irho = 1/rho;

xs2 = rho  .* ( xs*ct + ys*st );
ys2 = irho .* ( -xs*st + ys*ct );

r2 = xs2.^2 + ys2.^2;
g  = exp(-r2 / (2*sigmaSqr));

switch order
    case 0
        filter = g;
    case 1
        filter = -(xs2./sigmaSqr) .* g;
    case 2
        filter = (xs2.^2./sigmaSqr^2 - 1/sigmaSqr) .* g;
end
