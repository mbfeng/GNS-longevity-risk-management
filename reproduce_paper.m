%% REPRODUCE_PAPER  Reproduce all figures and tables of the paper.
%   Feng, Li & Zhou (2022), "Green nested simulation via likelihood ratio:
%   Applications to longevity risk management", Insurance: Mathematics and
%   Economics 106, 285-301.
%
%   RERUN_SIMULATIONS = false (default, about a minute): redraws Figs. 1-10 from the
%   saved results in results/, saves them to figures/, and prints Tables 1-3 with
%   Monte Carlo standard errors (also saved as results/imse_*.csv), then compares
%   them with the published values (verify_against_paper).
%   RERUN_SIMULATIONS = true: first regenerates results/ (about an hour on a 14-core
%   PC with the Parallel Computing Toolbox; without it the loops run serially).

RERUN_SIMULATIONS = false;

run(fullfile(fileparts(mfilename('fullpath')), 'setup_paths.m'));

if ~RERUN_SIMULATIONS && ~isfile(fullfile(fileparts(mfilename('fullpath')), 'results', 'case3_annuity.mat'))
    error(['results/*.mat not found. Download results.zip from the GitHub Release v1.0-ime2022 ' ...
           'and unzip it into results/, or set RERUN_SIMULATIONS = true (about an hour).']);
end

if RERUN_SIMULATIONS
    if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
        parpool('Threads');   % thread-based pool: light and needs no setup
    end
    case1_kcall_simulate(1000);
    case1_kcall_simulate(100);
    case2_qcallspread_simulate();
    case3_annuity_simulate();
end

fit_lee_carter();          % Fig. 1
case1_kcall_plots(1000);   % Figs. 2-4, Table 1
case1_kcall_plots(100);    % Fig. 5
case2_qcallspread_plots(); % Figs. 6, 7, 8a, Table 2
case3_annuity_plots();     % Figs. 8b, 9, 10, Table 3
verify_against_paper();    % Tables 1-3 vs the published values, within Monte Carlo error
