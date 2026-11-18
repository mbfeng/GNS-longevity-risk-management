function [I, se] = imse_table(Val, valTrue)
% IMSE_TABLE  Integrated mean squared errors, eq. (28), for every method and budget.
%   IMSE = mean over macro replications r of MSE_r, where MSE_r is the mean over
%   outer scenarios of (estimate - true value)^2 in replication r.
%
% Inputs:  Val (from run_methods), valTrue (M-by-1 true values)
% Outputs: I  : 4-by-B matrix; rows Nested, Approx, GNS, GNS_SN; columns = budgets
%          se : 4-by-B Monte Carlo standard errors of I (std of MSE_r / sqrt(R))

methods = {'Nested', 'Approx', 'GNS', 'GNS_SN'};
B = size(Val.Nested, 3);
R = size(Val.Nested, 2);
I = zeros(numel(methods), B);
se = zeros(numel(methods), B);
for j = 1:numel(methods)
    for b = 1:B
        mseR = mean((Val.(methods{j})(:, :, b) - valTrue).^2, 1);   % 1-by-R
        I(j, b) = mean(mseR);
        se(j, b) = std(mseR) / sqrt(R);
    end
end
end
