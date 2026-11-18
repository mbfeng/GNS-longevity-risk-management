function ax = plot_imse(I, budgets, labels, styles, ylims, figSize)
% PLOT_IMSE  IMSE versus simulation budget on a log scale (Figs. 4, 5 and 8).
%
% Inputs:
%   I       : methods-by-budgets matrix of IMSEs (rows of imse_table)
%   budgets : 1-by-B simulation budgets
%   labels  : legend label per method;  styles: line style per method, e.g. 'ko-'
%   ylims   : [ymin ymax];  figSize: [w h] in inches
% Output:
%   ax : axes handle (the paper's figures are exported from the axes only)

fig = figure('defaultTextInterpreter', 'latex');
B = numel(budgets);
args = {};
for j = 1:size(I, 1)
    args = [args, {1:B, I(j, :), styles{j}}]; %#ok<AGROW>
end
plt0 = plot(args{:}, 'LineWidth', 3, 'MarkerSize', 8);
ax = gca;
set(ax, 'XScale', 'linear', 'YScale', 'log', 'fontsize', 16)
axis([1, B, ylims])
legend(plt0, labels, 'location', 'southwest', 'fontsize', 16);
xticks(1:B)
xticklabels(arrayfun(@(b) sprintf('10^%d', round(log10(b))), budgets, 'UniformOutput', false))
ylabel('Integrated MSE', 'fontsize', 16)
xlabel('Simulation Budget', 'fontsize', 16)
pbaspect([1 1 1])
set(fig, 'Units', 'Inches', 'Position', [0.5, 0.5, figSize])
end
