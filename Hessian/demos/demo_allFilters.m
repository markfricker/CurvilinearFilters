% demo_allFilters.m
% Compare all Hessian-based filters on the shared reticulate test image.
% Saves composite figure to manual/figures/.

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..','src','engine'));
addpath(fullfile(thisDir,'..','src','utils'));
addpath(fullfile(thisDir,'..','tests','utilities'));
addpath(fullfile(thisDir,'..','..','Helper functions'));
figDir = fullfile(thisDir,'..','..','manual','figures');

I = im2single(makeReticulateTestImage(256, 256, 42));

filters = {'vesselness','ridge','blob','plate'};
titles  = {'Vesselness','Ridge','Blobness','Plateness'};

figure('Name','Hessian Filter Comparison','Color','w','Position',[100 100 1200 500]);

subplot(2,3,1);
imshow(I,[]);
title('Input');

for k = 1:numel(filters)
    [R,~,dir] = hessian2DFilters(I, ...
        'FilterType',filters{k});

    J = detectJunctions(R, dir);

    subplot(2,3,k+1);
    imshow(R,[]);
    hold on;
    [y,x] = find(J);
    plot(x,y,'ro','MarkerSize',8,'LineWidth',1.5);
    title(titles{k});
end

set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit',fullfile(figDir,'hessian_allFilters.pdf'));
fprintf('Saved: hessian_allFilters.pdf\n');

