function mask = vesselPolarityMask(Lambda1, WhiteOnDark)
%VESSELPOLARITYMASK Returns logical mask of valid vessel pixels
%
%   mask = vesselPolarityMask(Lambda1, WhiteOnDark)

if WhiteOnDark
    mask = (Lambda1 < 0);
else
    mask = (Lambda1 > 0);
end
end
