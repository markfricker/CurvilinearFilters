% demo_generateManualFigures.m
% Master script: generates all figures required by CurvilinearFilters.tex.
%
% Run this script from MATLAB before compiling the manual.
% All figures are saved to manual/figures/ as PDF files.
%
% Usage:
%   cd manual/demos
%   run demo_generateManualFigures

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
root    = fullfile(thisDir,'..');
figDir  = fullfile(root,'figures');

fprintf('=== CurvilinearFilters: generating manual figures ===\n\n');
fprintf('Output directory: %s\n\n', figDir);

% ---------------------------------------------------------------------------
% Shared paths
% ---------------------------------------------------------------------------
helperPath   = fullfile(root,'..','Helper functions');
hessianSrc   = fullfile(root,'..','Hessian','src');
hessianEng   = fullfile(root,'..','Hessian','src','engine');
hessianUtils = fullfile(root,'..','Hessian','src','utils');
hessianNeur  = fullfile(root,'..','Hessian','src','neuriteness');
mfatDrivers  = fullfile(root,'..','MFAT','src','drivers');
mfatCore     = fullfile(root,'..','MFAT','src','core');
mfatResp     = fullfile(root,'..','MFAT','src','responses');
mfatMod      = fullfile(root,'..','MFAT','src','modifiers');
mfatUtils    = fullfile(root,'..','MFAT','src','utils');
mfatCfg      = fullfile(root,'..','MFAT','config');
pcSrc        = fullfile(root,'..','PhaseCongruency','src');
pcUtils      = fullfile(root,'..','PhaseCongruency','utils');
soagkSrc     = fullfile(root,'..','SOAGK','src');
soagkUtils   = fullfile(root,'..','SOAGK','utils');
vessSrc      = fullfile(root,'..','Vesselness','src');
vessRoot     = fullfile(root,'..','Vesselness');

addpath(helperPath, hessianSrc, hessianEng, hessianUtils, hessianNeur, ...
        mfatDrivers, mfatCore, mfatResp, mfatMod, mfatUtils, mfatCfg, ...
        pcSrc, pcUtils, soagkSrc, soagkUtils, vessSrc, vessRoot);

% ---------------------------------------------------------------------------
% Shared test image
% ---------------------------------------------------------------------------
I = im2single(makeReticulateTestImage(256, 256, 42));

% ---------------------------------------------------------------------------
% 1. Hessian — composite (vesselness, ridge, blobness, plateness, neuriteness)
% ---------------------------------------------------------------------------
fprintf('[1/6] Hessian composite figure...\n');

[V,~,~] = hessian2DFilters(I,'FilterType','vesselness','Sigmas',1:6, ...
    'Parameters',struct('beta',0.5,'c',15));
[B,~,~] = hessian2DFilters(I,'FilterType','blob','Sigmas',1:6);
[P,~,~] = hessian2DFilters(I,'FilterType','plate','Sigmas',1:6, ...
    'Parameters',struct('alpha',0.5));
[Rd,~,~] = hessian2DFilters(I,'FilterType','ridge','Sigmas',1:6);
[N,~]   = neuriteness2D(I, 2);

fig = figure('Color','w','Position',[100 100 1400 500]);
tiledlayout(1,6,'Padding','compact','TileSpacing','compact');
tiles = {I, V, Rd, B, P, N};
tlbls = {'(a) Input','(b) Vesselness','(c) Ridge','(d) Blobness','(e) Plateness','(f) Neuriteness'};
for k = 1:6
    nexttile; imagesc(tiles{k}); axis image off;
    if k==1, colormap(gca,'gray'); else, colormap(gca,'hot'); end
    title(tlbls{k},'FontSize',8);
end
saveFig(fig, figDir, 'hessian_composite.pdf');

% ---------------------------------------------------------------------------
% 2. MFAT — lambda and prob side-by-side
% ---------------------------------------------------------------------------
fprintf('[2/6] MFAT figure...\n');

sigmas = 0.5:0.5:3;
R_lambda = mfatLambda(I, sigmas, 'tau',0.03,'tau2',0.3,'D',0.27, ...
    'whiteOnDark',true,'precision','single');
Prob = mfatProb(I, sigmas);

fig = figure('Color','w','Position',[100 100 900 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(I); axis image off; colormap(gca,'gray'); title('(a) Input','FontSize',8);
nexttile; imagesc(R_lambda); axis image off; colormap(gca,'hot'); title('(b) MFAT-\lambda','FontSize',8);
nexttile; imagesc(Prob,[0 1]); axis image off; colormap(gca,'hot'); title('(c) MFAT-Prob','FontSize',8);
saveFig(fig, figDir, 'mfat_result.pdf');

% ---------------------------------------------------------------------------
% 3. Phase Congruency
% ---------------------------------------------------------------------------
fprintf('[3/6] Phase congruency figure...\n');

[M, m, or_pc] = phaseCong3Optimised(I,'nscale',4,'norient',6, ...
    'minWaveLength',5,'mult',2.1,'sigmaOnf',0.55,'k',2.0);

fig = figure('Color','w','Position',[100 100 900 300]);
tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(I);     axis image off; colormap(gca,'gray'); title('(a) Input','FontSize',8);
nexttile; imagesc(M);     axis image off; colormap(gca,'hot');  title('(b) M (edges)','FontSize',8);
nexttile; imagesc(m);     axis image off; colormap(gca,'hot');  title('(c) m (corners)','FontSize',8);
nexttile; imagesc(or_pc); axis image off; colormap(gca,'hsv');  title('(d) Orientation','FontSize',8);
saveFig(fig, figDir, 'phaseCongruency_result.pdf');

% ---------------------------------------------------------------------------
% 4. Phase Symmetry
% ---------------------------------------------------------------------------
fprintf('[4/6] Phase symmetry figure...\n');

[phaseSym, orientation] = phaseSymOptimised(I,'nscale',5,'norient',6, ...
    'minWaveLength',3,'mult',2.1,'sigmaOnf',0.55,'k',2.0,'polarity',0);

fig = figure('Color','w','Position',[100 100 700 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(I);           axis image off; colormap(gca,'gray'); title('(a) Input','FontSize',8);
nexttile; imagesc(phaseSym);    axis image off; colormap(gca,'hot');  title('(b) Phase Symmetry','FontSize',8);
nexttile; imagesc(orientation); axis image off; colormap(gca,'hsv');  title('(c) Orientation','FontSize',8);
saveFig(fig, figDir, 'phaseSymmetry_result.pdf');

% ---------------------------------------------------------------------------
% 5. SOAGK
% ---------------------------------------------------------------------------
fprintf('[5/6] SOAGK figure...\n');

thetas = 0 : pi/6 : pi - pi/6;
[lineMap, dirMap] = AGlineDetectorSteerableConv2(I, 1:3, thetas, [1 2], 'single');

fig = figure('Color','w','Position',[100 100 700 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(I);       axis image off; colormap(gca,'gray'); title('(a) Input','FontSize',8);
nexttile; imagesc(lineMap); axis image off; colormap(gca,'hot');  title('(b) Line strength','FontSize',8);
nexttile; imagesc(dirMap);  axis image off; colormap(gca,'hsv');  title('(c) Orientation','FontSize',8);
saveFig(fig, figDir, 'soagk_result.pdf');

% ---------------------------------------------------------------------------
% 6. Vesselness (optimised standalone)
% ---------------------------------------------------------------------------
fprintf('[6/6] Optimised vesselness figure...\n');

opts_v.FrangiScaleRange = [1 6];
opts_v.FrangiScaleRatio = 2;
opts_v.FrangiBetaOne    = 0.5;
opts_v.FrangiBetaTwo    = 15;
opts_v.WhiteOnDark      = true;
opts_v.Precision        = 'single';
opts_v.verbose          = false;
[outIm, whatScale] = vesselnessFilter2DOptimised(I, opts_v);

fig = figure('Color','w','Position',[100 100 700 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(I);              axis image off; colormap(gca,'gray'); title('(a) Input','FontSize',8);
nexttile; imagesc(outIm);          axis image off; colormap(gca,'hot');  title('(b) Vesselness','FontSize',8);
nexttile; imagesc(double(whatScale)); axis image off; colormap(gca,'jet'); title('(c) Scale','FontSize',8);
saveFig(fig, figDir, 'vesselness_optimised_result.pdf');

% ---------------------------------------------------------------------------
fprintf('\nAll figures generated successfully.\n');
fprintf('Compile the manual with:\n');
fprintf('  cd manual && pdflatex CurvilinearFilters.tex\n');

% ---------------------------------------------------------------------------
function saveFig(fig, figDir, filename)
set(fig,'PaperPositionMode','auto');
outPath = fullfile(figDir, filename);
print(fig,'-dpdf','-bestfit',outPath);
fprintf('  Saved: %s\n', filename);
close(fig);
end
