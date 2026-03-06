% DEMO_MFATPROB  Probabilistic MFAT demonstration
% Uses the shared reticulate test image. Saves result to manual/figures/.

clear; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src','drivers'));
addpath(fullfile(thisDir,'..','src','core'));
addpath(fullfile(thisDir,'..','src','responses'));
addpath(fullfile(thisDir,'..','src','utils'));
addpath(fullfile(thisDir,'..','config'));
addpath(fullfile(thisDir,'..','..','Helper functions'));
figDir = fullfile(thisDir,'..','..','manual','figures');

% --- shared test image ---
I = im2single(makeReticulateTestImage(256, 256, 42));

sigmas = 0.5:0.5:3;

% --- MFAT-λ ---
R_lambda = mfatLambda(I, sigmas);

% --- MFAT-Prob ---
P = mfatProb(I, sigmas);

% --- visualize ---
figure('Color','w');
subplot(1,3,1);
imshow(I,[]); title('Input'); axis image off;

subplot(1,3,2);
imshow(R_lambda,[]); title('MFAT-\lambda'); axis image off;

subplot(1,3,3);
imshow(P,[]); title('MFAT-Prob (Posterior)'); axis image off;

set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'mfat_prob.pdf'));
fprintf('Saved: mfat_prob.pdf\n');
