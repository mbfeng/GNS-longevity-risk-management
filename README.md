# Code for "Green nested simulation via likelihood ratio: Applications to longevity risk management"

B.M. Feng, J.S.-H. Li and K.Q. Zhou, Insurance: Mathematics and Economics 106 (2022) 285-301,
https://doi.org/10.1016/j.insmatheco.2022.07.004. If you use this code, please cite the paper
(see `CITATION.cff`).
MATLAB code that reproduces every figure and table of the paper, within Monte Carlo error.

## Quick start

```matlab
reproduce_paper   % redraws Figs. 1-10 from results/, prints Tables 1-3 and checks them
                  % against the published values (about a minute)
```

The simulation results (`results/*.mat`, about 260 MB) are not in the git repository. Download
`results.zip` from the GitHub Release `v1.0-ime2022` and unzip it into `results/`, or regenerate them
by setting `RERUN_SIMULATIONS = true` in `reproduce_paper.m` (about an hour on a 14-core PC). Individual steps can also be run on their own, e.g. `case2_qcallspread_simulate()`
then `case2_qcallspread_plots()`.

Requirements: MATLAB (tested with R2023b) and the Statistics and Machine Learning Toolbox.
The Parallel Computing Toolbox is optional; without it the simulations run serially (slower; same results up to rounding, see below).

## Folder structure

```
Codes_Matlab/
  reproduce_paper.m   entry point
  setup_paths.m       adds functions/, scripts/ and external/plot_ci to the path
  scripts/            one function per case study and step (see the table below)
  functions/          building blocks (simulation, payoffs, estimators, plots)
  data/               hmd_ew_female_1962_2011.mat: Human Mortality Database deaths and
                      exposures, England & Wales females, ages 60-89, years 1962-2011
  results/            simulation results (*.mat) and IMSE tables with standard errors (*.csv)
  figures/            exported figures (PNG, named as in the manuscript)
  external/plot_ci/   third-party confidence-band plotting function (BSD license)
```

## Paper result -> code

| Paper result | Simulation | Figures and table |
|---|---|---|
| Fig. 1, kappa parameters | - | `fit_lee_carter` |
| Figs. 2-4, Table 1 (K-call, M = 1000) | `case1_kcall_simulate(1000)` | `case1_kcall_plots(1000)` |
| Fig. 5 (K-call, M = 100) | `case1_kcall_simulate(100)` | `case1_kcall_plots(100)` |
| Figs. 6, 7, 8a, Table 2 (q-call-spread) | `case2_qcallspread_simulate` | `case2_qcallspread_plots` |
| Figs. 8b, 9, 10, Table 3 (annuity) | `case3_annuity_simulate` | `case3_annuity_plots` |
| Tables 1-3 vs published values | - | `verify_against_paper` |

## How the simulations work

| Function | Role |
|---|---|
| `kappa_model` | Fits Lee-Carter (`lc_fit_svd`, `lc_fit_mle`) and the random walk with drift for kappa |
| `outer_scenarios` | M outer scenarios: equally spaced quantiles of kappa_tau |
| `simulate_kappa_paths` | Inner paths of kappa |
| `kcall_sim`, `qcallspread_sim`, `annuity_sim` | Discounted payoffs (and pathwise deltas) for each path |
| `kcall_exact`, `qcallspread_exact` | Exact values: closed form (K-call), numerical integration (q-call-spread) |
| `true_values_mc` | Monte Carlo "true" values for the annuity (independent paths per scenario) |
| `linear_taylor_sim`, `annuity_probit_taylor` | Taylor and probit-Taylor approximations |
| `gns_estimates` | GNS and self-normalized GNS with the mixture likelihood ratio, eq. (26) |
| `run_methods` | All methods x budgets x macro replications (parallel over replications) |
| `imse_table`, `write_imse_table` | IMSE, eq. (28), with Monte Carlo standard errors |
| `plot_setting`, `plot_bands`, `plot_imse` | The three kinds of figures |

Reproducibility: replication m always uses substream m of an `mrg32k3a` random stream with seed 2020,
so every run draws the same random numbers, with or without parallel workers. Rerunning on the same
machine and setup reproduces `results/` bit for bit. Serial and parallel runs agree to rounding
(GNS estimates differ by about 1e-13 relative, because multithreaded matrix products add in a
different order). Standard nested simulation and GNS use
the same inner paths in each replication (GNS reuses the nested-simulation paths, as in the paper).

## Code notation -> paper notation

| Code | Paper |
|---|---|
| `lc.a`, `lc.b`, `lc.k` | alpha_x, beta_x, kappa_t (ages 60-89 have indices 1-30) |
| `lc.kmu`, `lc.kVar` | drift and variance of the random walk for kappa |
| `lc.k0` | kappa_0 (2011) |
| `T`, `tau` | maturity and risk horizon, in years |
| `kOuter`, `M` | outer scenarios kappa_tau^(i) and their number |
| `khat` | anchor scenario of the (probit-)Taylor approximation |
| `NIn` = budget / M | N, inner paths per outer scenario in standard nested simulation |
| `budgets` | simulation budgets 1e3, 1e4, 1e5 |
| `NMac` | R, macro replications |
| `Val.Nested`, `Val.Approx`, `Val.GNS`, `Val.GNS_SN` | estimates of V by the four methods (M x R x 3) |

## Verification against the paper (2026-11-18, MATLAB R2023b)

The results in `results/` were regenerated with this code (about 50 minutes on 14 cores) and compared
with the 30 published entries of Tables 1-3 by `verify_against_paper`, which allows for the Monte
Carlo standard errors and for the paper's rounding to 3 digits. All 30 entries agree within 2.2
standard errors. Further checks:
- A clean copy of this folder, run with `RERUN_SIMULATIONS = true` in a fresh MATLAB session,
  reproduced the Case 1 and Case 2 results bit for bit. (Case 3 was regenerated afterwards, when the
  annuity true values were switched back to independent paths; see Notes.)
- A full rerun of Cases 1 and 2 with seed 2718 agrees with the paper within 2.4 standard errors, and
  the annuity true values computed with two other seeds give the same Table 3 to 3 digits.
- The regenerated figures match the published ones. Fig. 3a confirms that Table 1 and Figs. 2-4 use
  M = 1000 (1, 10, 100 inner simulations).
- An independent review of the code against the paper's equations and settings found no
  discrepancy affecting reported numbers.

`verify_against_paper` repeats this comparison for any set of results (it is the last step of
`reproduce_paper`). The original 2020-2022 results files were retired after this verification.

## Notes
- Section 4 of the paper says the "true" values use 10^7 inner paths; the published results used 10^6.
  Here the q-call-spread true values are computed exactly by numerical integration, and the annuity
  true values use 10^6 paths (`NTrue`), as in the published results.
- The annuity's "true" values use an independent set of 10^6 paths for every outer scenario, as in
  the published results. Do not replace this with one shared set of paths (common random numbers):
  the errors of the true values would then be one common shift instead of 1000 independent errors.
  That hardly matters for unbiased methods, but the probit-Taylor approximation has a bias of about
  0.004 of the same sign in almost every scenario, so its IMSE (Table 3, 99% bias at budget 1e5)
  would move by about 3%, several printed standard errors. (The printed standard errors cover only the
  R macro replications, not the error of the true values.) An intermediate version of this code did
  that; it gave 1.57e-5 or 1.64e-5 at budget 1e5 depending on the seed, against 1.62e-5 with
  independent paths and 1.63e-5 in the paper.
- Section 4.2 of the paper gives the q-call-spread's attachment and exhaustion points as
  AP = 3.703e-2 and EP = 5.037e-2. That exponent is a typo: the 5% and 95% quantiles of q_60,T are
  3.703e-3 and 5.037e-3, which is what the code and the published 2021 results use.
- The fitted kappa parameters match the paper: drift -0.4971, variance 0.7467 (volatility 0.8641),
  kappa_0 -14.2328.

## License
MIT (see `LICENSE`), for the code written for this project. `external/plot_ci/` keeps its own BSD
license. The data in `data/` come from the Human Mortality Database (University of California,
Berkeley, USA, and Max Planck Institute for Demographic Research, Germany; www.mortality.org)
and are subject to its terms of use.

## Change log
2026-11-19
- Added the MIT license and `CITATION.cff`; repository made public.

2026-11-18
- Restructured: scripts became functions (no `clearvars`), one entry point `reproduce_paper.m`, shared
  code moved into `functions/`, consistent snake_case names, unused options and outputs removed.
- Faster (about 1 hour instead of many hours for all results): parallel macro replications with
  reproducible random substreams, vectorized nested simulation, GNS reusing the nested-simulation
  paths, likelihood ratios computed once per block, exact q-call-spread true values by integration.
- Lee-Carter parameters are fitted from the HMD data on every run (no stored parameter file);
  the Econometrics Toolbox is no longer needed.
- Tables are printed with Monte Carlo standard errors and saved as CSV; legend labels are computed
  from the settings (the old Case 1 script showed the M = 100 labels).
- `setup_paths` no longer clears a caller's variable named `root`.
- Annuity true values: back to independent paths per scenario (as in the published results), after
  an intermediate shared-path version shifted the probit-Taylor IMSE by about 3%; Case 3 regenerated.
- Added `verify_against_paper` (compares Tables 1-3 with the published values, with z-scores).
- Verified (clean-copy rerun, second seed, figure comparison, independent code review against the
  paper), then deleted `results/original_runs/` (the 2020-2022 results, about 480 MB).

2026-11-17
- Folder cleaned up: removed an old copy of the Econometrics Toolbox, 2019 prototypes, the q-forward and
  y-option case studies (not in the paper) and a duplicate code folder.
