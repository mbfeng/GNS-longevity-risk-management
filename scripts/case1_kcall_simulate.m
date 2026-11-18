function res = case1_kcall_simulate(M, opts)
% CASE1_KCALL_SIMULATE  Case study 1: K-call (Section 4.1; Figs. 2-5, Table 1).
%   Values a K-call (maturity T = 10, strike K = E[kappa_T]) at the risk horizon
%   tau = 5 for M outer scenarios with standard nested simulation, Taylor
%   approximation, GNS and self-normalized GNS, each repeated NMac times.
%   True values are exact (closed form).
%
%   case1_kcall_simulate(1000)  Figs. 2-4 and Table 1 -> results/case1_kcall_M1000.mat
%   case1_kcall_simulate(100)   Fig. 5                 -> results/case1_kcall_M100.mat
%   case1_kcall_simulate(M, opts) overrides settings, e.g. for a quick test
%     opts = struct('NMac', 10, 'outFile', 'test.mat')
%
%   Draw the figures and the table with case1_kcall_plots(M).

if nargin < 1 || isempty(M), M = 1000; end
if nargin < 2, opts = struct(); end
p = project_paths();

% Settings (Section 4)
s.M = M;                                     % outer scenarios
s.seed = get_option(opts, 'seed', 2020);     % random seed
s.NMac = get_option(opts, 'NMac', 1000);     % macro replications R, eq. (28)
s.budgets = [1e3 1e4 1e5];                   % simulation budgets
s.NIn = s.budgets / M;                       % inner paths per outer scenario (nested sim.)
s.r = 0.03;                                  % discount rate, continuously compounded
s.T = 10;                                    % maturity, years
s.tau = 5;                                   % risk horizon, years
outFile = get_option(opts, 'outFile', fullfile(p.results, sprintf('case1_kcall_M%d.mat', M)));

lc = kappa_model();
s.K = lc.k0 + s.T*lc.kmu;                    % strike = E[kappa_T | kappa_0]
[kOuter, khat, outerPdf] = outer_scenarios(lc, s.tau, M);

% True values (closed form) and exact Taylor approximation (Fig. 2)
valTrue = kcall_exact(kOuter, lc.kmu, lc.kVol, s.K, s.T, s.tau, s.r);
[v0, d0] = kcall_exact(khat, lc.kmu, lc.kVol, s.K, s.T, s.tau, s.r);
valLin = v0 + d0 * (kOuter - khat);

% The four methods (Figs. 3-5, Table 1)
valueFcn = @(kProj) kcall_sim(kProj, s.K, s.T, s.tau, s.r);
approxFcn = @(kProj) linear_taylor_sim(kProj, kOuter, khat, valueFcn);
Val = run_methods(kOuter, khat, lc.kmu, lc.kVol, s.T - s.tau, s.budgets, s.NMac, s.seed, ...
    valueFcn, approxFcn);

res = struct('settings', s, 'lc', lc, 'kOuter', kOuter, 'khat', khat, 'outerPdf', outerPdf, ...
    'ValTrue', valTrue, 'ValLin', valLin, 'Val', Val);
save_results(outFile, res);
end
