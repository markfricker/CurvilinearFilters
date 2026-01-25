function [Gxx, Gyy, Gxy, pad] = getSOAGKSteerableBasisCached( ...
            sigma, rho, precision)

persistent cache

if isempty(cache)
    cache = containers.Map('KeyType','char','ValueType','any');
end

% --- robust key ---
key = sprintf('s=%.8g|r=%.8g|p=%s', sigma, rho, precision);

% --- cache hit ---
if isKey(cache, key)
    K   = cache(key);
    Gxx = K.Gxx;
    Gyy = K.Gyy;
    Gxy = K.Gxy;
    pad = K.pad;
    return
end

% --- cache miss ---
params.rho  = rho;
params.size = round(9 * sigma * max(1, 1/rho));
if mod(params.size,2)==0
    params.size = params.size + 1;
end

[Gxx, Gyy, Gxy] = createSOAGKSteerableBasisConv2( ...
                    sigma, params, precision);

pad = floor(size(Gxx,1)/2);

% --- store ---
K.Gxx = Gxx;
K.Gyy = Gyy;
K.Gxy = Gxy;
K.pad = pad;

cache(key) = K;
end
