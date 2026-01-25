% demo_junctionRecovery.m
% Demonstrate junction recovery on a network

clear; clc;

I = generateJunctionTestImage(256);

[R,~,dir] = hessian2DFilters(I, ...
    'FilterType','vesselness');

J = detectJunctions(R, dir, ...
    'Radius',3, ...
    'ThresholdFrac',0.2);

figure;
tiledlayout
nexttile
imshow(I,[]);
title('Input');

nexttile
imshow(R,[]);
imshow(R,[]);
hold on;
[y,x] = find(J);
plot(x,y,'ro','MarkerSize',8,'LineWidth',1.5);
title('Recovered Junction Nodes');
