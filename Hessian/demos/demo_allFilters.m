% demo_allFilters.m
% Compare all Hessian-based filters

clear; clc;

I = generateTestImage(256);

filters = {'vesselness','ridge','blob','plate'};
titles  = {'Vesselness','Ridge','Blobness','Plateness'};

figure('Name','Hessian Filter Comparison','Position',[100 100 1200 500]);

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

