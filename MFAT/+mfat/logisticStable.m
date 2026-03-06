function p = logisticStable(x, realClass)
% logisticStable  numerically-stable logistic function, result in realClass
    if nargin < 2, realClass = 'single'; end
    xmax = cast(20, realClass);
    xmin = cast(-20, realClass);
    xc = min(max(x, xmin), xmax);
    p = cast(1, realClass) ./ (cast(1, realClass) + exp(-xc));
end
