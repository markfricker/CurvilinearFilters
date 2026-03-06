% demo_basicVesselness.m
% Basic vesselness example using the shared reticulate test image.
% Saves result to manual/figures/.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..','src','engine'));
addpath(fullfile(thisDir,'..','..','Helper functions'));
figDir = fullfile(thisDir,'..','..','manual','figures');

I = im2single(makeReticulateTestImage(256, 256, 42));

[R, scale, dir] = hessian2DFilters(I, ...
    'FilterType','vesselness', ...
    'WhiteOnDark',true);

figure('Color','w');
subplot(1,3,1);
imshow(I,[]);
title('Input');

subplot(1,3,2);
imshow(R,[]);
title('Vesselness Response');

subplot(1,3,3);
imshow(dir,[]);
title('Orientation (rad)');

set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'hessian_vesselness.pdf'));
fprintf('Saved: hessian_vesselness.pdf\n');
