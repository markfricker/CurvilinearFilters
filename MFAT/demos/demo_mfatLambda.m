% DEMO_MFATLAMBDA  Demonstration of deterministic MFAT-λ
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

% --- MFAT scales ---
sigmas = 0.5:0.5:3;

% --- run MFAT-λ ---
R_lambda = mfatLambda(I, sigmas, ...
    'tau',0.03, ...
    'tau2',0.3, ...
    'D',0.27, ...
    'whiteOnDark',true, ...
    'precision','single');

% --- display ---
figure('Color','w');
subplot(1,2,1);
imshow(I,[]); title('Input'); axis image off;
subplot(1,2,2);
imshow(R_lambda,[]); title('MFAT-\lambda (Deterministic)'); axis image off;

set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'mfat_lambda.pdf'));
fprintf('Saved: mfat_lambda.pdf\n');
