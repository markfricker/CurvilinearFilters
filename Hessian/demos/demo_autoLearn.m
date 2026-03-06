% demo_autoLearn.m
% Demonstrate automatic scale selection on the shared reticulate test image.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..','src','engine'));
addpath(fullfile(thisDir,'..','src','utils'));
addpath(fullfile(thisDir,'..','..','Helper functions'));

I = im2single(makeReticulateTestImage(256, 256, 42));

filters = {'vesselness','ridge','blob','plate'};
titles  = {'Vesselness','Ridge','Blobness','Plateness'};

sigmas = autoLearnHessianScales(I, 0.5, 8, 0.5);

figure('Name','Hessian Filter Comparison (Auto-Scales)','Position',[100 100 1200 500]);

subplot(2,3,1);
imshow(I,[]);
title(['Input - sigmas ' num2str(sigmas)]);

for k = 1:numel(filters)
    R = hessian2DFilters(I, ...
        'FilterType',filters{k},'Sigmas',sigmas);

    % Recompute eigenvalues at a representative scale
    [Dxx,Dxy,Dyy] = applyHessian2D(single(I), 2);
    [l1,l2,~,~]   = eig2image(Dxx,Dxy,Dyy);
    C = hessianConfidence(R,l1,l2);
    subplot(2,3,k+1);
    imshowpair(R,C,'montage');
    title(titles{k});
end
