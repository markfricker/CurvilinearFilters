function test_vesselnessFilter2DOptimised()
% Unit test for FrangiFilter2dOptimised
%
% Validates numerical equivalence against the reference
% FrangiFilter2D implementation within tolerance.

fprintf('Running FrangiFilter2dOptimised unit test...\n');

% -------------------------------------------------------------
% Synthetic test image (tubular structure)
% -------------------------------------------------------------
N = 256;
[x,y] = meshgrid(1:N,1:N);
cx = N/2; cy = N/2;

radius = 3;
I = exp(-((x-cx).^2 + (y-cy).^2)/(2*radius^2));
I = I + 0.05*randn(size(I));   % add noise
I = mat2gray(I);

% -------------------------------------------------------------
% Options
% -------------------------------------------------------------
options.FrangiScaleRange = [1 6];
options.FrangiScaleRatio = 1;
options.FrangiBetaOne = 0.5;
options.FrangiBetaTwo = 15;
options.WhiteOnDark = true;
options.Precision = 'double';
options.verbose = false;

% -------------------------------------------------------------
% Run reference implementation
% -------------------------------------------------------------
tic
[outRef, scaleRef, dirRef] = FrangiFilter2D(double(I), options);
toc
% -------------------------------------------------------------
% Run optimised implementation
% -------------------------------------------------------------
tic
[outOpt, scaleOpt, dirOpt] = vesselnessFilter2DOptimised(double(I), options);
toc

figure
tiledlayout(2,3)
nexttile
imshow(outRef,[])
title 'frangi response'
nexttile
imshow(scaleRef,[])
title 'frangi scale'
nexttile
imshow(dirRef,[])
title 'frangi dir'
nexttile
imshow(outOpt,[])
title 'opt response'
nexttile
imshow(scaleOpt,[])
title 'optscale'
nexttile
imshow(dirOpt,[])
title 'opt dir'

figure
% -------------------------------------------------------------
% Assertions
% -------------------------------------------------------------
tolResponse = 1e-10;
tolAngle    = 1e-10;

assert(max(abs(outRef(:) - outOpt(:))) < tolResponse, ...
    'Vesselness response mismatch');
imshowpair(scaleRef,scaleOpt)

assert(isequal(scaleRef, scaleOpt), ...
    'Scale selection mismatch');

% Direction only valid where vesselness > 0
mask = outRef > 0;
assert(max(abs(dirRef(mask) - dirOpt(mask))) < tolAngle, ...
    'Direction mismatch');

fprintf('✔ FrangiFilter2dOptimised PASSED all tests\n');
end
