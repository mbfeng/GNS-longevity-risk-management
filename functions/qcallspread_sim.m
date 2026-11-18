function [value, delta] = qcallspread_sim(kProj, a, b, x, T, tau, r, AP, EP)
% QCALLSPREAD_SIM  Discounted payoffs and pathwise deltas of a q-call-spread (Section 2.3.2).
%   Payoff at T: min(max((q - AP)/(EP - AP), 0), 1), where
%   q = 1 - exp(-exp(a_x + b_x*kappa_T)) is the death probability at the reference age.
%
% Inputs:
%   kProj : (T-tau)-by-N matrix of kappa paths after tau (last row = kappa_T)
%   a, b  : Lee-Carter alpha and beta, indexed by age (1 = age 60)
%   x     : index of the reference age
%   T     : maturity, in years from time 0
%   tau   : risk horizon, in years
%   r     : discount rate, continuously compounded, per year
%   AP    : attachment point
%   EP    : exhaustion point
% Outputs:
%   value : N-by-1 discounted payoffs
%   delta : N-by-1 pathwise derivatives in kappa_tau

T2M = T-tau;
ax = a(x);
bx = b(x);
disc = exp(-r*T2M);
yT = (ax+bx*kProj(T2M,:))';     % log death rate at T
eY = exp(yT);
eMinusW = exp(-eY);
qx = 1 - eMinusW;               % death probability at T

A = (qx - AP)./(EP - AP);       % payoff before capping at 0 and 1
A(A<0) = 0;
A(A>1) = 1;
value = disc.*A;

Ind = (0<A) & (A<1);
qx_delta = bx.*eY.*eMinusW;     % derivative of q in kappa_T
delta = disc.*(qx_delta.*Ind)'./(EP - AP);
end
