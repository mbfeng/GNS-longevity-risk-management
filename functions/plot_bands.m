function fig = plot_bands(kOuter, valTrue, V, labels, o)
% PLOT_BANDS  95% bands of an estimator over the macro replications, one band per
% simulation budget, with the true value (panels of Figs. 3, 7 and 10).
%
% Inputs:
%   kOuter  : M-by-1 outer scenarios
%   valTrue : M-by-1 true values
%   V       : M-by-nMac-by-B estimates (one field of the output of run_methods)
%   labels  : 1-by-B legend labels of the budgets
%   o       : options struct: yrng, axPos, legendLoc, figSize [w h] in inches,
%             trueLineWidth

colors = [186 228 179; 116 196 118; 35 139 69] / 255;   % light to dark green
xrng = [min(kOuter), max(kOuter)];
fig = figure('defaultTextInterpreter', 'latex');
hold on
B = size(V, 3);
patches = gobjects(B, 1);
for b = 1:B
    Y = [mean(V(:, :, b), 2), quantile(V(:, :, b), [0.025, 0.975], 2)];
    h = plot_ci(kOuter, Y, 'PatchColor', colors(b, :), 'PatchAlpha', 0.75, ...
        'MainLineColor', 'none', 'LineColor', 'none');
    patches(b) = h.Patch;
end
plt0 = plot(kOuter, valTrue, 'LineWidth', o.trueLineWidth, 'Color', 'k');
legend([plt0; patches], [{'True Value'}, labels], 'location', o.legendLoc, 'fontsize', 20);
xlabel('Outer kappa values, $\kappa_{\tau}$', 'fontsize', 20)
axis([xrng, o.yrng])
set(gca, 'fontsize', 20, 'position', o.axPos, 'box', 'on');
hold off
set(fig, 'Units', 'Inches', 'Position', [0.5, 0.5, o.figSize])
end
