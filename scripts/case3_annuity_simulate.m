function res = case3_annuity_simulate(opts)
% CASE3_ANNUITY_SIMULATE  Case study 3: temporary life annuity (Section 4.3; Figs. 8b, 9, 10, Table 3).
%   A 30-year temporary life annuity issued at age 60 is valued at the risk horizon
%   tau = 10 (age 70, 20 years left) for M = 1000 outer scenarios with standard nested
%   simulation, the probit-Taylor approximation and (self-normalized) GNS, each
%   repeated NMac times. True values: Monte Carlo with NTrue paths.
%
%   case3_annuity_simulate()      -> results/case3_annuity.mat
%   case3_annuity_simulate(opts)  overrides settings, e.g. for a quick test
%     opts = struct('NMac', 10, 'NTrue', 1e4, 'M', 100, 'outFile', 'test.mat')
%
%   Draw the figures and the table with case3_annuity_plots().

if nargin < 1, opts = struct(); end
p = project_paths();

% Settings (Section 4)
s.M = get_option(opts, 'M', 1000);           % outer scenarios
s.seed = get_option(opts, 'seed', 2020);     % random seed
s.NMac = get_option(opts, 'NMac', 1000);     % macro replications R, eq. (28)
s.NTrue = get_option(opts, 'NTrue', 1e6);    % paths for the "true" values. Section 4 of
                                             % the paper says 10^7; the published Table 3
                                             % used 1e6.
s.budgets = [1e3 1e4 1e5];                   % simulation budgets
s.NIn = s.budgets / s.M;                     % inner paths per outer scenario (nested sim.)
s.r = 0.03;                                  % discount rate, continuously compounded
s.T = 30;                                    % term of the annuity, years
s.tau = 10;                                  % risk horizon, years
s.xTau = 1 + s.tau;                          % index of the age at tau: 11 = age 70
outFile = get_option(opts, 'outFile', fullfile(p.results, 'case3_annuity.mat'));

lc = kappa_model();
[kOuter, khat, outerPdf] = outer_scenarios(lc, s.tau, s.M);
nSteps = s.T - s.tau;
valueFcn = @(kProj) annuity_sim(kProj, lc.a, lc.b, s.xTau, s.T, s.tau, s.r);
approxFcn = @(kProj) annuity_probit_taylor(kProj, khat, kOuter, lc.a, lc.b, s.xTau, s.T, s.tau, s.r);

% "True" values and probit-Taylor approximation from NTrue paths (Fig. 9)
[valTrue, valLin] = true_values_mc(kOuter, khat, lc.kmu, lc.kVol, nSteps, s.NTrue, s.seed + 1, ...
    valueFcn, approxFcn);

% The methods (Figs. 8b, 10, Table 3)
Val = run_methods(kOuter, khat, lc.kmu, lc.kVol, nSteps, s.budgets, s.NMac, s.seed, ...
    valueFcn, approxFcn);

res = struct('settings', s, 'lc', lc, 'kOuter', kOuter, 'khat', khat, 'outerPdf', outerPdf, ...
    'ValTrue', valTrue, 'ValLin', valLin, 'Val', Val);
save_results(outFile, res);
end
