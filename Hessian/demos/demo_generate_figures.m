% demo_generate_figures.m
% Generates all filter response figures from the complex phantom

addpath('../core','../engine','../neuriteness','../tests/utils');

% -------------------------------------------------------------------------
% Load phantom (PNG, not PDF)
% -------------------------------------------------------------------------
I = im2double(imread('../manual/figures/complex_phantom.png'));

% -------------------------------------------------------------------------
% Vesselness
% -------------------------------------------------------------------------
[V,~,~] = hessian2DFilters(I, ...
    'FilterType','vesselness', ...
    'Sigmas',1:6, ...
    'Parameters',struct('beta',0.5,'c',15));

figure; imagesc(V); axis image off; colormap hot;
title('Vesselness response');
print('-dpdf','../manual/figures/complex_vesselness.pdf');

% -------------------------------------------------------------------------
% Blobness
% -------------------------------------------------------------------------
[B,~,~] = hessian2DFilters(I, ...
    'FilterType','blob', ...
    'Sigmas',1:6);

figure; imagesc(B); axis image off; colormap hot;
title('Blobness response');
print('-dpdf','../manual/figures/complex_blobness.pdf');

% -------------------------------------------------------------------------
% Plateness
% -------------------------------------------------------------------------
[P,~,~] = hessian2DFilters(I, ...
    'FilterType','plate', ...
    'Sigmas',1:6, ...
    'Parameters',struct('alpha',0.5));

figure; imagesc(P); axis image off; colormap hot;
title('Plateness response');
print('-dpdf','../manual/figures/complex_plateness.pdf');

% -------------------------------------------------------------------------
% Neuriteness (single scale)
% -------------------------------------------------------------------------
[N,~] = neuriteness2D(I,2);

figure; imagesc(N); axis image off; colormap hot;
title('Neuriteness response');
print('-dpdf','../manual/figures/complex_neuriteness.pdf');
