function demo_generate_composite_figure()
%DEMO_GENERATE_COMPOSITE_FIGURE
% Generates a composite figure showing:
%  (a) Input phantom
%  (b) Vesselness
%  (c) Blobness
%  (d) Plateness
%  (e) Neuriteness
%
% Output: manual/figures/complex_composite.pdf

addpath('../core','../engine','../neuriteness','../tests/utils');

% -------------------------------------------------------------------------
% Load phantom
% -------------------------------------------------------------------------
I = im2double(imread('../manual/figures/complex_phantom.png'));

% -------------------------------------------------------------------------
% Compute responses
% -------------------------------------------------------------------------
[V,~,~] = hessian2DFilters(I, ...
    'FilterType','vesselness', ...
    'Sigmas',1:6, ...
    'Parameters',struct('beta',0.5,'c',15));

[B,~,~] = hessian2DFilters(I, ...
    'FilterType','blob', ...
    'Sigmas',1:6);

[P,~,~] = hessian2DFilters(I, ...
    'FilterType','plate', ...
    'Sigmas',1:6, ...
    'Parameters',struct('alpha',0.5));

[N,~] = neuriteness2D(I,2);

% -------------------------------------------------------------------------
% Composite figure
% -------------------------------------------------------------------------
figure('Color','w','Position',[100 100 1200 800]);

tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% (a) Input
nexttile;
imagesc(I); axis image off; colormap gray;
title('(a) Input phantom');

% (b) Vesselness
nexttile;
imagesc(V); axis image off; colormap hot;
title('(b) Vesselness');

% (c) Blobness
nexttile;
imagesc(B); axis image off; colormap hot;
title('(c) Blobness');

% (d) Plateness
nexttile;
imagesc(P); axis image off; colormap hot;
title('(d) Plateness');

% (e) Neuriteness
nexttile;
imagesc(N); axis image off; colormap hot;
title('(e) Neuriteness');

% Empty tile (intentional, for balance)
nexttile;
axis off;

% -------------------------------------------------------------------------
% Save composite figure (fit to page)
% -------------------------------------------------------------------------
set(gcf,'PaperPositionMode','auto');
print('-dpdf','-bestfit','../manual/figures/complex_composite.pdf');
close;