function [value, delta] = kcall_sim(kProj, K, T, tau, r)
% KCALL_SIM  Discounted payoffs and pathwise deltas of a K-call from simulated kappa paths.
%
% Inputs:
%   kProj : (T-tau)-by-N matrix of kappa paths after tau (last row = kappa_T)
%   K     : strike on kappa
%   T     : maturity, in years from time 0
%   tau   : risk horizon, in years
%   r     : discount rate, continuously compounded, per year
% Outputs:
%   value : N-by-1 discounted payoffs
%   delta : N-by-1 pathwise derivatives in kappa_tau

T2M = T-tau;
disc = exp(-r*T2M);

A = kProj(T2M, :) - K;
Ind = A>0;
value = disc.*(A.*Ind)';
delta = disc.*Ind';
end
