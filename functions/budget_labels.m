function labels = budget_labels(budgets, nIn)
% BUDGET_LABELS  Legend labels for the simulation budgets, e.g. '10^3 sim. budget'.
%   With nIn (inner paths per outer scenario) the labels also give the number of
%   inner simulations, e.g. '1 inner sim, 10^3 sim. budget' (Fig. 3a).
if nargin < 2
    labels = arrayfun(@(b) sprintf('10^%d sim. budget', round(log10(b))), ...
        budgets, 'UniformOutput', false);
else
    labels = arrayfun(@(n, b) sprintf('%d inner sim, 10^%d sim. budget', n, round(log10(b))), ...
        nIn, budgets, 'UniformOutput', false);
end
end
