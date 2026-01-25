% demo_generateDocumentationFigures.m
% Generate figures for LaTeX documentation

clear; clc;

outDir = fullfile('doc','figures');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% --- Generate test image ---
I = generateTestImage(256);
imwrite(I, fullfile(outDir,'test_image.png'));

% --- Vesselness ---
V = hessian2DFilters(I, 'FilterType','vesselness');
imwrite(mat2gray(V), fullfile(outDir,'vesselness.png'));

% --- Ridge ---
R = hessian2DFilters(I, 'FilterType','ridge');
imwrite(mat2gray(R), fullfile(outDir,'ridge.png'));

% --- Neuriteness ---
N = hessian2DFilters(I, 'FilterType','neuriteness');
imwrite(mat2gray(N), fullfile(outDir,'neuriteness.png'));

% --- Blobness ---
B = hessian2DFilters(I, 'FilterType','blob');
imwrite(mat2gray(B), fullfile(outDir,'blob.png'));

% --- Plateness ---
P = hessian2DFilters(I, 'FilterType','plate');
imwrite(mat2gray(P), fullfile(outDir,'plate.png'));

fprintf('Documentation figures written to %s\n', outDir);
