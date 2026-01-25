% DEMO_MFATLAMBDA  Demonstration of deterministic MFAT-λ

clear; close all;

% --- load example image ---
I = im2single(imread('example.png'));   % replace with your image
if size(I,3) > 1
    I = rgb2gray(I);
end

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
figure;
imshow(R_lambda,[]);
title('MFAT-\lambda (Deterministic)');
