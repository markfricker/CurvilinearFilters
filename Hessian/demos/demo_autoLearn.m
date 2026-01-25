% demo_allFilters.m
% Compare all Hessian-based filters

clear; clc;

I = generateTestImage(256);

filters = {'vesselness','ridge','neuriteness','blob','plate'};
titles  = {'Vesselness','Ridge','Neuriteness','Blobness','Plateness'};


sigmas = autoLearnHessianScales(I,0.5,8,0.5);

figure('Name','Hessian Filter Comparison','Position',[100 100 1200 500]);

subplot(2,3,1);
imshow(I,[]);
title(['Input ' num2str(sigmas)]);

for k = 1:numel(filters)
    R = hessian2DFilters(I, ...
        'FilterType',filters{k},'Sigmas',sigmas);

    % Recompute eigenvalues at a representative scale
    [Dxx,Dxy,Dyy] = Hessian2D(single(I), 2);
    [l1,l2,~,~]   = eig2image(Dxx,Dxy,Dyy);
    C = hessianConfidence(R,l1,l2);
    subplot(2,3,k+1);
    imshowpair(R,C,'montage');
    title(titles{k});
end

