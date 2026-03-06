% demo_phaseSymmetry.m
% Demonstrates phaseSymOptimised on the shared reticulate test image.
% Saves symmetry/orientation figure to manual/figures/.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..','utils'));
addpath(fullfile(thisDir,'..','..','Helper functions'));
figDir = fullfile(thisDir,'..','..','manual','figures');

% --- shared test image ---
I = im2single(makeReticulateTestImage(256, 256, 42));

% --- phase symmetry ---
[phaseSym, orientation, totalEnergy] = phaseSymOptimised(I, ...
    'nscale',     5, ...
    'norient',    6, ...
    'minWaveLength', 3, ...
    'mult',       2.1, ...
    'sigmaOnf',   0.55, ...
    'k',          2.0, ...
    'polarity',   0);

% --- display ---
figure('Name','Phase Symmetry','Color','w','Position',[100 100 900 300]);

subplot(1,4,1);
imshow(I,[]); title('Input'); axis image off;

subplot(1,4,2);
imshow(phaseSym,[]); title('Phase Symmetry'); axis image off; colormap(gca,'hot');

subplot(1,4,3);
imagesc(orientation); axis image off; colormap(gca,'hsv');
title('Orientation'); colorbar;

subplot(1,4,4);
imshow(totalEnergy,[]); title('Total Energy'); axis image off; colormap(gca,'hot');

% --- save ---
set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'phaseSymmetry_result.pdf'));
fprintf('Saved: phaseSymmetry_result.pdf\n');
