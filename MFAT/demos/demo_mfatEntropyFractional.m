% DEMO_MFATENTROPYFRACTIONAL  Optional MFAT modifiers
% Uses the shared reticulate test image.

clear; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src','drivers'));
addpath(fullfile(thisDir,'..','src','core'));
addpath(fullfile(thisDir,'..','src','responses'));
addpath(fullfile(thisDir,'..','src','modifiers'));
addpath(fullfile(thisDir,'..','src','utils'));
addpath(fullfile(thisDir,'..','config'));
addpath(fullfile(thisDir,'..','..','Helper functions'));

% --- shared test image ---
I = im2single(makeReticulateTestImage(256, 256, 42));

sigmas = 0.5:0.5:3;

% --- MFAT-λ ---
R_lambda = mfatLambda(I, sigmas);

% --- collect FA stack ---
opts = mfatConfig().core;
faStack = zeros([size(I), numel(sigmas)], 'single');

for k = 1:numel(sigmas)
    geom = mfatCore2D(I, sigmas(k), opts);
    faStack(:,:,k) = geom.fa;
end

% --- entropy weighting ---
R_ent = mfatEntropyWeight(R_lambda, faStack, ...
    'beta',0.5, ...
    'combine','mean');

% --- fractional shaping ---
R_frac = mfatFractional(R_lambda, 'alpha',0.7);

% --- display ---
figure('Color','w');
subplot(1,3,1); imshow(R_lambda,[]); title('MFAT-\lambda');
subplot(1,3,2); imshow(R_ent,[]);    title('Entropy-weighted');
subplot(1,3,3); imshow(R_frac,[]);   title('Fractional (\alpha=0.7)');
