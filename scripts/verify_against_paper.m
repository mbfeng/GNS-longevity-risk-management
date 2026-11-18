function T = verify_against_paper(resultsDir)
% VERIFY_AGAINST_PAPER  Compare the IMSEs in results/ with the published Tables 1-3.
%   verify_against_paper()     checks results/ (e.g. after RERUN_SIMULATIONS = true)
%   verify_against_paper(dir)  checks the results files in another folder
%
%   For every published entry, prints the regenerated IMSE, the paper's value and
%   z = (regenerated - paper) / sqrt(2*se^2 + h^2/3), where se is the Monte Carlo standard
%   error of the regenerated IMSE over the R macro replications (the paper's value is an
%   estimate with about the same standard error, hence 2*se^2) and h is half a unit in the
%   last of the 3 digits printed in the paper (rounding error, uniform on +-h). |z| below
%   about 2.5 means agreement within Monte Carlo and rounding error.
%
%   se does not include the Monte Carlo error of the annuity's "true" values. That error
%   matters only for the biased probit-Taylor row of Table 3, and only slightly, because
%   the true values use independent paths per scenario (see README.md, Notes).
%
%   Returns a table with one row per published entry.

p = project_paths();
if nargin < 1 || isempty(resultsDir), resultsDir = p.results; end

% Published values (columns: budgets 1e3, 1e4, 1e5). Rows: nested, Taylor, GNS, GNS-SN;
% NaN where the paper does not report the method.
paper.case1_kcall_M1000 = [1.08e0  1.08e-1 1.08e-2     % Table 1
                           1.84e-1 1.83e-1 1.83e-1
                           1.55e-2 1.55e-3 1.49e-4
                           4.47e-3 4.92e-4 4.58e-5];
paper.case2_qcallspread = [2.90e-2 2.91e-3 2.90e-4     % Table 2
                           4.01e-4 3.69e-4 3.66e-4
                           NaN     NaN     NaN
                           7.35e-5 7.27e-6 7.63e-7];
paper.case3_annuity     = [2.30e-2 2.30e-3 2.30e-4     % Table 3
                           3.71e-5 1.81e-5 1.63e-5
                           NaN     NaN     NaN
                           1.27e-4 1.23e-5 1.28e-6];
methodNames = {'Nested', 'Taylor', 'GNS', 'GNS_SN'};
budgetNames = {'1e3', '1e4', '1e5'};

[caseCol, methodCol, budgetCol] = deal({});
[regen, pub, se, z] = deal([]);
files = fieldnames(paper);
for c = 1:numel(files)
    S = load(fullfile(resultsDir, [files{c} '.mat']), 'Val', 'ValTrue');
    [I, s] = imse_table(S.Val, S.ValTrue);
    P = paper.(files{c});
    for j = 1:4
        for b = 1:3
            if isnan(P(j, b)), continue, end
            caseCol{end+1, 1} = files{c};          %#ok<AGROW>
            methodCol{end+1, 1} = methodNames{j};  %#ok<AGROW>
            budgetCol{end+1, 1} = budgetNames{b};  %#ok<AGROW>
            regen(end+1, 1) = I(j, b);             %#ok<AGROW>
            pub(end+1, 1) = P(j, b);               %#ok<AGROW>
            se(end+1, 1) = s(j, b);                %#ok<AGROW>
            % the paper prints 3 significant digits: rounding error is uniform on +-h
            h = 0.5 * 10^(floor(log10(P(j, b))) - 2);
            z(end+1, 1) = (I(j, b) - P(j, b)) / sqrt(2 * s(j, b)^2 + h^2 / 3);   %#ok<AGROW>
        end
    end
end
T = table(caseCol, methodCol, budgetCol, regen, pub, se, z, 'VariableNames', ...
    {'case', 'method', 'budget', 'regenerated', 'paper', 'se', 'z'});

fprintf('\nRegenerated IMSEs vs the published Tables 1-3 (%s)\n', resultsDir);
fprintf('%-20s %-7s %-6s %12s %12s %8s\n', 'case', 'method', 'budget', 'regenerated', 'paper', 'z');
for i = 1:height(T)
    flag = '';
    if abs(T.z(i)) > 2.5, flag = '  <-- |z| > 2.5'; end
    fprintf('%-20s %-7s %-6s %12.3g %12.3g %8.1f%s\n', T.case{i}, T.method{i}, T.budget{i}, ...
        T.regenerated(i), T.paper(i), T.z(i), flag);
end
nOut = sum(abs(T.z) > 2.5);
fprintf('%d of %d entries within 2.5 standard errors of the paper (largest |z| = %.1f).\n', ...
    height(T) - nOut, height(T), max(abs(T.z)));
end
