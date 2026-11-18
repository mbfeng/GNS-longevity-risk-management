function I = case1_kcall_plots(M, export)
% CASE1_KCALL_PLOTS  Case study 1 (K-call): figures and IMSE table.
%   case1_kcall_plots(1000)  Fig. 2, Figs. 3a-d, Fig. 4 and Table 1
%   case1_kcall_plots(100)   Fig. 5 (and the same figures for M = 100)
%   case1_kcall_plots(M, false) draws the figures without saving them.
%
%   Reads results/case1_kcall_M<M>.mat; saves figures to figures/ and the IMSE table
%   to results/imse_case1_kcall_M<M>.csv. Returns the 4-by-3 IMSE matrix.

if nargin < 1 || isempty(M), M = 1000; end
if nargin < 2, export = true; end
p = project_paths();
S = load(fullfile(p.results, sprintf('case1_kcall_M%d.mat', M)));
s = S.settings;
sfx = '';
if M ~= 1000, sfx = sprintf('_M%d', M); end

% Fig. 2: problem setting
f = plot_setting(S.kOuter, S.ValTrue, S.ValLin, S.outerPdf, struct('yrng', [-2 6], ...
    'axPos', [0.05 0.15 0.9 0.8], 'approxLabel', 'Taylor approx.', 'legendLoc', 'northwest'));
if export, export_figure(f, ['KOptionSetting' sfx]); end

% Figs. 3a-d: 95% bands
o = struct('yrng', [-2 8], 'axPos', [0.06 0.11 0.9 0.85], 'legendLoc', 'northwest', ...
    'figSize', [10 7], 'trueLineWidth', 2.5);
panels = {'Nested', 'KOptionNestedSim', budget_labels(s.budgets, s.NIn)
          'Approx', 'KOptionApprox',    budget_labels(s.budgets)
          'GNS',    'KOptionGNS',       budget_labels(s.budgets)
          'GNS_SN', 'KOptionGNSSN',     budget_labels(s.budgets)};
for j = 1:size(panels, 1)
    f = plot_bands(S.kOuter, S.ValTrue, S.Val.(panels{j, 1}), panels{j, 3}, o);
    if export, export_figure(f, [panels{j, 2} sfx]); end
end

% Table 1 and Fig. 4 (M = 1000) or Fig. 5 (M = 100): IMSE
[I, se] = imse_table(S.Val, S.ValTrue);
names = {'Standard nested sim.', 'Taylor approx.', 'Green nested sim.', 'GNS, self-normalized'};
if M == 1000
    tableTitle = 'Table 1: IMSE, K-call, M = 1000';
else
    tableTitle = sprintf('Fig. 5 data: IMSE, K-call, M = %d', M);
end
write_imse_table(I, se, 1:4, names, s.budgets, tableTitle, ...
    fullfile(p.results, sprintf('imse_case1_kcall_M%d.csv', M)));
ax = plot_imse(I, s.budgets, {'Standard nested sim.', 'Taylor approximation', ...
    'Green nested sim.', 'GNS, self-normalized'}, {'ko-', 'rs:', 'm^-.', 'b^--'}, [1e-5 1e1], [10 8]);
if export, export_figure(ax, sprintf('KOptionIMSE%d', M)); end
end
