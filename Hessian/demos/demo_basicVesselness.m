% demo_basicVesselness.m
% Basic vesselness example

clear; clc;

I = im2single(imread('cameraman.tif'));

[R, scale, dir] = hessian2DFilters(I, ...
    'FilterType','vesselness', ...
    'WhiteOnDark',true);

figure;
subplot(1,3,1);
imshow(I,[]);
title('Input');

subplot(1,3,2);
imshow(R,[]);
title('Vesselness Response');

subplot(1,3,3);
imshow(dir,[]);
title('Orientation (rad)');
