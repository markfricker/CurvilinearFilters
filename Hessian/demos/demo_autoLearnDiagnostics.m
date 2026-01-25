% demo_autolearn_diagnostics.m
% Visualise how auto-learn selects scales

clear; clc;

I = generateTestImage(256);

sigmaMin  = 0.5;
sigmaMax  = 8;
sigmaStep = 0.5;

% --- Diagnostic plot ---
plotAutoLearnDiagnostics(I, sigmaMin, sigmaMax, sigmaStep);

% --- Use auto-learned sigmas ---
sigmas = autoLearnHessianScales(I, sigmaMin, sigmaMax, sigmaStep);

fprintf('Auto-learn selected sigmas:\n');
disp(sigmas);

% --- Apply a filter using these sigmas ---
R = hessian2DFilters(I, ...
    'FilterType','vesselness', ...
    'Sigmas',sigmas);

figure;
imshow(R,[]);
title('Vesselness with Auto-Learned Scales');
