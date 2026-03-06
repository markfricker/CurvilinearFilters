% demo_phaseCongruency.m
% Demonstrates phaseCong3Optimised on the shared reticulate test image.
% Saves edge/corner/orientation figure to manual/figures/.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..','utils'));
addpath(fullfile(thisDir,'..','..','Helper functions'));
figDir = fullfile(thisDir,'..','..','manual','figures');

% --- shared test image ---
I = im2single(makeReticulateTestImage(256, 256, 42));

% --- phase congruency ---
[M, m, or] = phaseCong3Optimised(I, ...
    'nscale',    4, ...
    'norient',   6, ...
    'minWaveLength', 5, ...
    'mult',      2.1, ...
    'sigmaOnf',  0.55, ...
    'k',         2.0);

% --- display ---
figure('Name','Phase Congruency','Color','w','Position',[100 100 900 300]);

subplot(1,4,1);
imshow(I,[]); title('Input'); axis image off;

subplot(1,4,2);
imshow(M,[]); title('M (edges)'); axis image off; colormap(gca,'hot');

subplot(1,4,3);
imshow(m,[]); title('m (corners)'); axis image off; colormap(gca,'hot');

subplot(1,4,4);
imagesc(or); axis image off; colormap(gca,'hsv');
title('Orientation'); colorbar;

% --- save ---
set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'phaseCongruency_result.pdf'));
fprintf('Saved: phaseCongruency_result.pdf\n');
