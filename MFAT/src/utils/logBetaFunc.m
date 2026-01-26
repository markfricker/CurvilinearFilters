function lb = logBetaFunc(a, b, realClass)
% logBetaFunc  compute log(Beta(a,b)) robustly, return in realClass
    lbDouble = gammaln(double(a)) + gammaln(double(b)) - gammaln(double(a + b));
    if nargin < 3, realClass = 'single'; end
    lb = cast(lbDouble, realClass);
end
