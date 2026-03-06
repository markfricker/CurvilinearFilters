function test_vesselnessFilter2DOptimised()
% Unit test for FrangiFilter2dOptimised
%
% Validates numerical equivalence against the reference
% FrangiFilter2D implementation within tolerance.

fprintf('Running FrangiFilter2dOptimised unit test...\n');

% -------------------------------------------------------------
% Synthetic test image (tubular structure)
% -------------------------------------------------------------
% N = 256;
% [x,y] = meshgrid(1:N,1:N);
% cx = N/2; cy = N/2;
% 
% radius = 3;
% I = exp(-((x-cx).^2 + (y-cy).^2)/(2*radius^2));
% I = I + 0.05*randn(size(I));   % add noise
% I = mat2gray(I);
I = makeReticulateTestImage(256, 256, 123); % deterministic synthetic

% -------------------------------------------------------------
% Options
% -------------------------------------------------------------
options.FrangiScaleRange = [1 6];
options.FrangiScaleRatio = 1;
options.FrangiBetaOne = 0.5;
options.FrangiBetaTwo = 15;
options.BlackWhite = false;    % FrangiFilter2D field: false = white ridges
options.WhiteOnDark = true;    % vesselnessFilter2DOptimised field (same polarity)
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

% -------------------------------------------------------------
% Assertions
% -------------------------------------------------------------
tolResponse = 1e-10;
tolAngle    = 1e-10;

assert(max(abs(outRef(:) - outOpt(:))) < tolResponse, ...
    'Vesselness response mismatch');
imshowpair(scaleRef,scaleOpt)

% Scale maps must agree at nearly all pixels; allow up to 2% mismatch
% (tie-breaking at equal-response scales can differ by one step).
fracMismatch = mean(scaleRef(:) ~= scaleOpt(:));
assert(fracMismatch < 0.02, ...
    sprintf('Scale selection mismatch: %.1f%% of pixels differ', 100*fracMismatch));

% Direction only valid where vesselness > 0
mask = outRef > 0;
assert(max(abs(dirRef(mask) - dirOpt(mask))) < tolAngle, ...
    'Direction mismatch');

fprintf('✔ FrangiFilter2dOptimised PASSED all tests\n');
end
