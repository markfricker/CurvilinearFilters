function mfatDisplayCompare(I, outLam, outProb, outFrac, outEnt)
% Robust MFAT comparison display (no saturation bias)

figure('Color','w','Name','MFAT Comparison (Corrected)');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% ---- input ----
nexttile;
imagesc(I); axis image off;
colormap gray;
title('Input');

% robust normalization for non-probabilistic outputs
normResp = @(x) ...
    max(0, min(1, ...
    (x - prctile(x(:),5)) / ...
    (prctile(x(:),95) - prctile(x(:),5) + eps)));

% ---- MFAT-λ ----
nexttile;
imagesc(normResp(outLam));
axis image off;
colormap gray;
title('MFAT-\lambda');

% ---- MFAT-Prob ----
nexttile;
imagesc(outProb,[0 1]);
axis image off;
colormap gray;
title('MFAT-Prob');

% ---- MFAT-Frac ----
nexttile;
imagesc(normResp(outFrac));
axis image off;
colormap gray;
title('MFAT-Frac');

% ---- MFAT-Entropy ----
nexttile;
imagesc(normResp(outEnt));
axis image off;
colormap gray;
title('MFAT-Entropy');

% ---- colorbar ----
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb,'Normalized response');
end
