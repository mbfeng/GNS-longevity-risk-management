function export_figure(h, name)
% EXPORT_FIGURE  Save a figure or axes as figures/<name>.png at 300 dpi.
p = project_paths();
exportgraphics(h, fullfile(p.figures, [name '.png']), 'Resolution', 300);
end
