function out = imcomplementSafe(x, realClass)
% imcomplementSafe - class-preserving imcomplement (1 - x) cast-safe
    if nargin < 2, realClass = 'single'; end
    out = cast(1, realClass) - cast(x, realClass);
end
