function res = case2_qcallspread_simulate(opts)
% CASE2_QCALLSPREAD_SIMULATE  Case study 2: q-call-spread (Section 4.2; Figs. 6, 7, 8a, Table 2).
%   Values a q-call-spread on the death probability q at age 60 at maturity T = 10,
%   at the risk horizon tau = 5, for M = 1000 outer scenarios with standard nested
%   simulation, Taylor approximation and (self-normalized) GNS, each repeated NMac times.
%   True values and the Taylor approximation of Fig. 6 are computed by numerical
%   integration (the paper used Monte Carlo with many paths; the difference is
%   negligible).
%
%   case2_qcallspread_simulate()      -> results/case2_qcallspread.mat
%   case2_qcallspread_simulate(opts)  overrides settings, e.g. for a quick test
%     opts = struct('NMac', 10, 'M', 100, 'outFile', 'test.mat')
%
%   Draw the figures and the table with case2_qcallspread_plots().

if nargin < 1, opts = struct(); end
p = project_paths();

% Settings (Section 4)
s.M = get_option(opts, 'M', 1000);           % outer scenarios
s.seed = get_option(opts, 'seed', 2020);     % random seed
s.NMac = get_option(opts, 'NMac', 1000);     % macro replications R, eq. (28)
s.budgets = [1e3 1e4 1e5];                   % simulation budgets
s.NIn = s.budgets / s.M;                     % inner paths per outer scenario (nested sim.)
s.r = 0.03;                                  % discount rate, continuously compounded
s.T = 10;                                    % maturity, years
s.tau = 5;                                   % risk horizon, years
s.x = 1;                                     % index of the reference age: 1 = age 60
outFile = get_option(opts, 'outFile', fullfile(p.results, 'case2_qcallspread.mat'));

lc = kappa_model();

% Attachment and exhaustion points: 5% and 95% quantiles of q at T
kT = norminv([0.05 0.95], lc.k0 + s.T*lc.kmu, sqrt(s.T*lc.kVar));
qT = 1 - exp(-exp(lc.a(s.x) + lc.b(s.x)*kT));
s.AP = qT(1);
s.EP = qT(2);

[kOuter, khat, outerPdf] = outer_scenarios(lc, s.tau, s.M);

% True values and exact Taylor approximation (Fig. 6), by numerical integration
exactFcn = @(k) qcallspread_exact(k, lc.a, lc.b, s.x, s.T, s.tau, s.r, s.AP, s.EP, lc.kmu, lc.kVol);
valTrue = exactFcn(kOuter);
[v0, d0] = exactFcn(khat);
valLin = v0 + d0 * (kOuter - khat);

% The methods (Figs. 7, 8a, Table 2)
valueFcn = @(kProj) qcallspread_sim(kProj, lc.a, lc.b, s.x, s.T, s.tau, s.r, s.AP, s.EP);
approxFcn = @(kProj) linear_taylor_sim(kProj, kOuter, khat, valueFcn);
Val = run_methods(kOuter, khat, lc.kmu, lc.kVol, s.T - s.tau, s.budgets, s.NMac, s.seed, ...
    valueFcn, approxFcn);

res = struct('settings', s, 'lc', lc, 'kOuter', kOuter, 'khat', khat, 'outerPdf', outerPdf, ...
    'ValTrue', valTrue, 'ValLin', valLin, 'Val', Val);
save_results(outFile, res);
end
