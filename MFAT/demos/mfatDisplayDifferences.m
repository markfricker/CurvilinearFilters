function mfatDisplayDifferences(outLam, outProb, outFrac, outEnt)
% mfatDisplayDifferences  Difference maps vs MFAT-lambda

figure('Color','w','Name','MFAT Differences');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% symmetric color range
dMax = max(abs([ ...
    outProb(:)-outLam(:);
    outFrac(:)-outLam(:);
    outEnt(:)-outLam(:)]));

% ---- Prob - Lambda ----
nexttile;
imagesc(outProb - outLam,[-dMax dMax]);
axis image off;
colormap redblue;
title('MFAT-Prob − MFAT-\lambda');

% ---- Frac - Lambda ----
nexttile;
imagesc(outFrac - outLam,[-dMax dMax]);
axis image off;
colormap redblue;
title('MFAT-Frac − MFAT-\lambda');

% ---- Entropy - Lambda ----
nexttile;
imagesc(outEnt - outLam,[-dMax dMax]);
axis image off;
colormap redblue;
title('MFAT-Entropy − MFAT-\lambda');

colorbar;
end
