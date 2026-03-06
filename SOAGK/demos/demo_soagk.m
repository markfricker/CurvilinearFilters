% demo_soagk.m
% Demonstrates AGlineDetectorSteerableConv2 on the shared reticulate test image.
% Saves line strength and orientation figure to manual/figures/.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..','utils'));
addpath(fullfile(thisDir,'..','..','Helper functions'));
figDir = fullfile(thisDir,'..','..','manual','figures');

% --- shared test image ---
I = im2single(makeReticulateTestImage(256, 256, 42));

% --- SOAGK parameters ---
sigmas = [1 2 3];
thetas = 0 : pi/6 : pi - pi/6;   % 6 orientations in [0, pi)
rhos   = [1 2];                    % isotropic and elongated

% --- run detector ---
[lineMap, dirMap] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos, 'single');

% --- display ---
figure('Name','SOAGK Line Detector','Color','w','Position',[100 100 900 300]);

subplot(1,3,1);
imshow(I,[]); title('Input'); axis image off;

subplot(1,3,2);
imshow(lineMap,[]); title('Line strength'); axis image off; colormap(gca,'hot');

subplot(1,3,3);
imagesc(dirMap); axis image off; colormap(gca,'hsv');
title('Orientation (rad)'); colorbar;

% --- save ---
set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'soagk_result.pdf'));
fprintf('Saved: soagk_result.pdf\n');
