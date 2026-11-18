function I = case3_annuity_plots(export)
% CASE3_ANNUITY_PLOTS  Case study 3 (temporary life annuity): Fig. 9, Figs. 10a-c, Fig. 8b, Table 3.
%   case3_annuity_plots()       draws and saves the figures
%   case3_annuity_plots(false)  draws the figures without saving them
%
%   Reads results/case3_annuity.mat; saves figures to figures/ and the IMSE table
%   to results/imse_case3_annuity.csv. Returns the 4-by-3 IMSE matrix.

if nargin < 1, export = true; end
p = project_paths();
S = load(fullfile(p.results, 'case3_annuity.mat'));
s = S.settings;
yrng = [11.25, 13];

% Fig. 9: problem setting, with a zoomed inset at kappa = -11
inset = struct('center', [-11, 11.54], 'size', [0.45, 0.11], 'axPos', [0.75 0.6 0.2 0.3]);
f = plot_setting(S.kOuter, S.ValTrue, S.ValLin, S.outerPdf, struct('yrng', yrng, ...
    'axPos', [0.075 0.15 0.9 0.8], 'approxLabel', 'Linear approx.', 'legendLoc', 'west', ...
    'inset', inset));
if export, export_figure(f, 'KAnnuitySetting'); end

% Figs. 10a-c: 95% bands
o = struct('yrng', yrng, 'axPos', [0.1 0.1 0.85 0.85], 'legendLoc', 'southwest', ...
    'figSize', [6 9], 'trueLineWidth', 2.5);
panels = {'Nested', 'KAnnuityNestedSim', 2.5; 'Approx', 'KAnnuityApprox', 1.5
          'GNS_SN', 'KAnnuityGNSSN', 1.5};
for j = 1:size(panels, 1)
    o.trueLineWidth = panels{j, 3};
    f = plot_bands(S.kOuter, S.ValTrue, S.Val.(panels{j, 1}), budget_labels(s.budgets), o);
    if export, export_figure(f, panels{j, 2}); end
end

% Table 3 and Fig. 8b: IMSE
[I, se] = imse_table(S.Val, S.ValTrue);
rows = [1 2 4];
write_imse_table(I, se, rows, {'Standard nested sim.', '(Probit-)Taylor approx.', 'GNS, self-normalized'}, ...
    s.budgets, 'Table 3: IMSE, temporary life annuity', fullfile(p.results, 'imse_case3_annuity.csv'));
ax = plot_imse(I(rows, :), s.budgets, {'Standard nested sim.', 'Taylor approximation', ...
    'GNS, self-normalized'}, {'ko-', 'rs:', 'b^--'}, [5e-7 1e-1], [9 8]);
if export, export_figure(ax, 'KAnnuityIMSE'); end
end
