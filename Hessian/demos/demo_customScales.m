% demo_customScales.m
% Effect of sigma range on Hessian response

clear; clc;

I = generateTestImage(256);

sigmaSets = {
    0.5:0.5:2,   'Small scales';
    1:1:6,       'Medium scales';
    3:1:8,       'Large scales'
};

figure('Name','Scale Sensitivity','Position',[100 100 1200 400]);

for k = 1:size(sigmaSets,1)
    sigmas = sigmaSets{k,1};
    label  = sigmaSets{k,2};

    R = hessian2DFilters(I, ...
        'FilterType','vesselness', ...
        'Sigmas',sigmas);

    subplot(1,3,k);
    imshow(R,[]);
    title(label);
end
