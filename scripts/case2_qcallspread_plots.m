function I = case2_qcallspread_plots(export)
% CASE2_QCALLSPREAD_PLOTS  Case study 2 (q-call-spread): Fig. 6, Figs. 7a-c, Fig. 8a, Table 2.
%   case2_qcallspread_plots()       draws and saves the figures
%   case2_qcallspread_plots(false)  draws the figures without saving them
%
%   Reads results/case2_qcallspread.mat; saves figures to figures/ and the IMSE table
%   to results/imse_case2_qcallspread.csv. Returns the 4-by-3 IMSE matrix.

if nargin < 1, export = true; end
p = project_paths();
S = load(fullfile(p.results, 'case2_qcallspread.mat'));
s = S.settings;

% Fig. 6: problem setting
f = plot_setting(S.kOuter, S.ValTrue, S.ValLin, S.outerPdf, struct('yrng', [-0.2 1], ...
    'axPos', [0.1 0.15 0.85 0.8], 'approxLabel', 'Taylor approx.', 'legendLoc', 'northwest'));
if export, export_figure(f, 'qCallSpreadSetting'); end

% Figs. 7a-c: 95% bands (GNS without self-normalization is not in the paper)
o = struct('yrng', [-0.2 1], 'axPos', [0.1 0.1 0.85 0.85], 'legendLoc', 'northwest', ...
    'figSize', [6 9], 'trueLineWidth', 2.5);
panels = {'Nested', 'qCallSpreadNestedSim'; 'Approx', 'qCallSpreadApprox'
          'GNS', 'qCallSpreadGNS'; 'GNS_SN', 'qCallSpreadGNSSN'};
for j = 1:size(panels, 1)
    f = plot_bands(S.kOuter, S.ValTrue, S.Val.(panels{j, 1}), budget_labels(s.budgets), o);
    if export, export_figure(f, panels{j, 2}); end
end

% Table 2 and Fig. 8a: IMSE
[I, se] = imse_table(S.Val, S.ValTrue);
rows = [1 2 4];
write_imse_table(I, se, rows, {'Standard nested sim.', 'Taylor approx.', 'GNS, self-normalized'}, ...
    s.budgets, 'Table 2: IMSE, q-call-spread', fullfile(p.results, 'imse_case2_qcallspread.csv'));
ax = plot_imse(I(rows, :), s.budgets, {'Standard nested sim.', 'Taylor approximation', ...
    'GNS, self-normalized'}, {'ko-', 'rs:', 'b^--'}, [5e-7 1e-1], [9 8]);
if export, export_figure(ax, 'qCallSpreadIMSE'); end
end
