% DEMO_MFATPROB  Probabilistic MFAT demonstration

clear; close all;

I = im2single(imread('example.png'));
if size(I,3) > 1
    I = rgb2gray(I);
end

sigmas = 0.5:0.5:3;

% --- MFAT-λ ---
R_lambda = mfatLambda(I, sigmas);

% --- MFAT-Prob ---
P = mfatProb(I, sigmas);

% --- visualize ---
figure;
subplot(1,2,1);
imshow(R_lambda,[]);
title('MFAT-\lambda');

subplot(1,2,2);
imshow(P,[]);
title('MFAT-Prob (Posterior)');
