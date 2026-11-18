function value = annuity_probit_taylor(kProj, khat, kOuter, a, b, x, T, tau, r)
% ANNUITY_PROBIT_TAYLOR  Probit-Taylor approximation of a temporary life annuity (Lemma 2.3, eq. (17)).
%   Each survival probability S_t is approximated by Phi(c0 + c1*(kappa_tau - khat))
%   around the anchor scenario khat. The coefficients c0, c1 are estimated from the
%   paths in kProj, which start at khat.
%
% Inputs:
%   kProj  : (T-tau)-by-N matrix of kappa paths after tau, starting from khat
%   khat   : anchor scenario
%   kOuter : M-by-1 outer scenarios at which the annuity is valued
%   a, b   : Lee-Carter alpha and beta, indexed by age (1 = age 60)
%   x      : index of the annuitant's age at tau
%   T      : end of the annuity, in years from time 0
%   tau    : risk horizon, in years
%   r      : discount rate, continuously compounded, per year
% Output:
%   value  : M-by-1 approximate annuity values

T2M = T-tau;
disc = exp(-r*(1:T2M));        % discount factors

ax = a(x:(x+T2M-1));
bx = b(x:(x+T2M-1));

Y = ax+bx.*kProj;              % log death rates
eY = exp(Y);
W = cumsum(eY, 1);
eMinusW = exp(-W);             % survival probabilities on each path
temp1 = cumsum(bx.*eY, 1);

% Coefficients of the probit-Taylor expansion, one pair per annuity payment:
% c0 = probit of the survival probability, c1 = its derivative in kappa.
Lia0 = norminv(mean(eMinusW,2));
Del0 = -mean(eMinusW.*temp1, 2) ./ normpdf(Lia0);

L_Fvalue = Lia0 + Del0*(kOuter' - khat);
L_Pvalue = normcdf(L_Fvalue);  % approximate survival probabilities
value = (disc*L_Pvalue)';
end
