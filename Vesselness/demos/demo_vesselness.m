% demo_vesselness.m
% Demonstrates the optimised standalone Frangi vesselness filter on the
% shared reticulate test image.
% Saves result figure to manual/figures/.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..'));              % for +vesselness package
addpath(fullfile(thisDir,'..','..','Hessian','src'));    % applyHessian2D, eig2image
addpath(fullfile(thisDir,'..','..','Helper functions')); % makeReticulateTestImage
figDir = fullfile(thisDir,'..','..','manual','figures');

% --- shared test image ---
I = im2single(makeReticulateTestImage(256, 256, 42));

% --- vesselness options ---
options.FrangiScaleRange = [1 6];
options.FrangiScaleRatio = 2;
options.FrangiBetaOne    = 0.5;
options.FrangiBetaTwo    = 15;
options.WhiteOnDark      = true;
options.Precision        = 'single';
options.verbose          = false;

% --- run filter ---
[outIm, whatScale, Direction] = vesselnessFilter2DOptimised(I, options);

% --- display ---
figure('Name','Optimised Vesselness','Color','w','Position',[100 100 900 300]);

subplot(1,3,1);
imshow(I,[]); title('Input'); axis image off;

subplot(1,3,2);
imshow(outIm,[]); title('Vesselness'); axis image off; colormap(gca,'hot');

subplot(1,3,3);
imagesc(double(whatScale)); axis image off; colormap(gca,'jet');
title('Scale index'); colorbar;

% --- save ---
set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'vesselness_optimised_result.pdf'));
fprintf('Saved: vesselness_optimised_result.pdf\n');
