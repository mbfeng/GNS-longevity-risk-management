function [value, delta, gamma] = annuity_sim(kProj, a, b, x, T, tau, r)
% ANNUITY_SIM  Value of a temporary life annuity at tau for each simulated kappa path.
%   The annuity pays 1 at the end of each year from tau+1 to T if the annuitant is
%   alive. For each path, value = sum_t exp(-r*t) * S_t, where
%   S_t = exp(-sum_{s<=t} exp(a_{x+s-1} + b_{x+s-1}*kappa_{tau+s})) is the survival probability.
%
% Inputs:
%   kProj : (T-tau)-by-N matrix of kappa paths after tau (row t = year tau+t)
%   a, b  : Lee-Carter alpha and beta, indexed by age (1 = age 60)
%   x     : index of the annuitant's age at tau
%   T     : end of the annuity, in years from time 0
%   tau   : risk horizon, in years
%   r     : discount rate, continuously compounded, per year
% Outputs:
%   value : N-by-1 annuity values
%   delta, gamma : N-by-1 pathwise first and second derivatives in kappa_tau
%                  (only computed if requested; not used for the paper)

T2M = T-tau;
disc = exp(-r*(1:T2M));      % discount factors

Y = a(x:(x+T2M-1)) + b(x:(x+T2M-1)) .* kProj(1:T2M, :);
eY = exp(Y);
W = cumsum(eY, 1);

value = (disc*exp(-W))';
if nargout > 1
    delta = -(disc*(exp(-W).*cumsum(b(x:(x+T2M-1)).*eY, 1)))';
    gamma = (disc*(exp(-W).*(cumsum(b(x:(x+T2M-1)).*eY, 1).^2 - cumsum(b(x:(x+T2M-1)).^2.*eY, 1))))';
end
end
