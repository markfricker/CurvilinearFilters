function plotAutoLearnDiagnostics(I, sigmaMin, sigmaMax, sigmaStep)
% plotAutoLearnDiagnostics
% Diagnostic plot for Hessian auto-learn scale selection
%
% Shows:
%  1) Scale-normalised Hessian energy vs sigma
%  2) Relative energy change (ΔE)
%  3) Selected sigmas

I = single(I);
sigmas = sigmaMin:sigmaStep:sigmaMax;

energy = zeros(size(sigmas),'single');

% ---- Compute energy ----
for k = 1:numel(sigmas)
    s = sigmas(k);
    [Dxx,~,Dyy] = Hessian2D(I, s);
    curv = abs(s^2 * Dxx) + abs(s^2 * Dyy);
    energy(k) = mean(curv(:));
end

% ---- Normalise energy ----
energy = energy / max(energy);

% ---- Relative change ----
dE = [0 diff(energy) ./ energy(1:end-1)];

% ---- Auto-learn selection ----
keep = false(size(sigmas));
keep(1) = true;
keep(2:end) = dE(2:end) > 0.05;

% ---- Plot ----
figure('Name','Auto-learn Diagnostics','Position',[100 100 800 600]);

subplot(2,1,1);
plot(sigmas, energy, '-o','LineWidth',1.5);
hold on;
plot(sigmas(keep), energy(keep), 'ro','MarkerSize',8,'LineWidth',1.5);
xlabel('\sigma');
ylabel('Normalised Hessian Energy');
title('Scale-normalised Hessian Energy');
grid on;
legend('Energy','Selected scales','Location','SouthEast');

subplot(2,1,2);
plot(sigmas, dE, '-s','LineWidth',1.5);
hold on;
yline(0.05,'r--','Threshold (5%)','LineWidth',1.2);
xlabel('\sigma');
ylabel('Relative Energy Change (\DeltaE)');
title('Information Gain Across Scales');
grid on;

end
