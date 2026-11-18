function T = write_imse_table(I, se, rows, rowNames, budgets, title, csvFile)
% WRITE_IMSE_TABLE  Print an IMSE table in the layout of Tables 1-3 and save it as CSV.
%
% Inputs:
%   I, se    : IMSEs and their standard errors (from imse_table)
%   rows     : indices of the rows of I to report
%   rowNames : names of those rows, as in the paper
%   budgets  : simulation budgets (column headers)
%   title    : printed title
%   csvFile  : output file (full precision, with standard errors)
% Output:
%   T : the table

b = arrayfun(@(x) sprintf('1e%d', round(log10(x))), budgets, 'UniformOutput', false);
T = array2table([I(rows, :), se(rows, :)], 'RowNames', rowNames, ...
    'VariableNames', [strcat('imse_', b), strcat('se_', b)]);
writetable(T, csvFile, 'WriteRowNames', true);

fprintf('\n%s  (IMSE, Monte Carlo standard error in brackets)\n%-24s', title, 'Simulation budget');
fprintf('%22s', b{:});
fprintf('\n');
for j = 1:numel(rows)
    fprintf('%-24s', rowNames{j});
    fprintf('%11.3g (%7.2g)', [I(rows(j), :); se(rows(j), :)]);
    fprintf('\n');
end
end
