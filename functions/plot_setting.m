function fig = plot_setting(kOuter, valTrue, valApprox, outerPdf, o)
% PLOT_SETTING  "Problem settings" figure (Figs. 2, 6, 9): true value, Taylor
% approximation without simulation error, and the density of the outer scenarios.
%
% Options (struct o): yrng, axPos (axes position), approxLabel, legendLoc,
%   inset (optional struct with fields center [x y], size [w h], axPos) for a zoomed
%   inset as in Fig. 9.

xrng = [min(kOuter), max(kOuter)];
fig = figure('defaultTextInterpreter', 'latex');
axes('position', o.axPos, 'YColor', 'k', 'box', 'off', 'fontsize', 20);
yyaxis left
p1 = line(kOuter, valTrue, 'LineWidth', 3, 'LineStyle', '-', 'Color', '#0072BD');
p2 = line(kOuter, valApprox, 'LineWidth', 3, 'LineStyle', ':', 'Color', '#D95319');
axis([xrng, o.yrng])
hasInset = isfield(o, 'inset') && ~isempty(o.inset);
if hasInset   % mark the zoomed region
    c = o.inset.center; s = o.inset.size;
    rectangle('Position', [c - s/2, s], 'EdgeColor', [0.4, 0.1, 0.4], 'LineWidth', 1);
end
yyaxis right
a1 = area(kOuter, outerPdf, 0, 'EdgeColor', 'none', 'FaceColor', [0.4660, 0.6740, 0.1880], 'FaceAlpha', .5);
ylim([0 0.7])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'none';
ax.Box = 'off';
legend([p1; p2; a1], {'True value', o.approxLabel, 'Outer kappa pdf'}, 'location', o.legendLoc, 'fontsize', 20)
xlabel('Outer kappa values, $\kappa_{\tau}$', 'fontsize', 20)
if hasInset   % zoomed inset
    axes('Units', 'Normalized', 'Position', o.inset.axPos, 'XTick', [], 'YTick', [], ...
        'Box', 'on', 'LineWidth', 2, 'Color', [1, 1, 1]);
    hold on
    plot(kOuter, valTrue, '-', kOuter, valApprox, ':', 'LineWidth', 3)
    xlabel(sprintf('detail at $\\kappa = %g$', c(1)), 'fontsize', 20);
    axis([c(1) - s(1)/2, c(1) + s(1)/2, c(2) - s(2)/2, c(2) + s(2)/2]);
    hold off
end
set(fig, 'Units', 'Inches', 'Position', [0.5, 0.5, 10, 5])
end
