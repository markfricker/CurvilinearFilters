function I = generateTestBlobImage(sz, radius)
% generateTestBlobImage  Deterministic sz×sz image with a centred Gaussian blob.
%
%   I = generateTestBlobImage(sz, radius)
%
%   sz     - image side length in pixels (default 64)
%   radius - blob radius in pixels; sigma = radius/2 (default 5)
%
%   Returns a single-precision [0,1] image.  No toolbox required.

if nargin < 1, sz     = 64; end
if nargin < 2, radius = 5;  end

[x, y] = meshgrid(1:sz, 1:sz);
cx  = (sz + 1) / 2;
cy  = (sz + 1) / 2;
sig = radius / 2;
I   = single(exp(-((x - cx).^2 + (y - cy).^2) / (2 * sig^2)));
I   = I / max(I(:));
end
