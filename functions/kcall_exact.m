function [value, delta] = kcall_exact(kTau, kmu, kVol, K, T, tau, r)
% KCALL_EXACT  Closed-form value and delta at tau of a K-call (Section 2.4.1).
%   The K-call pays max(kappa_T - K, 0) at T; kappa follows the random walk with drift.
%
% Inputs:
%   kTau : scalar or vector of kappa_tau values (outer scenarios)
%   kmu  : drift of kappa per year
%   kVol : volatility of kappa per year
%   K    : strike on kappa
%   T    : maturity, in years from time 0
%   tau  : risk horizon, in years
%   r    : discount rate, continuously compounded, per year
% Outputs (same size as kTau):
%   value : value at tau
%   delta : derivative of the value in kappa_tau

T2M = T-tau;                % remaining time to maturity
disc = exp(-r*T2M);         % discount factor
kVolT2M = kVol*sqrt(T2M);   % volatility over the remaining time to maturity

A = K - (kTau + T2M*kmu);
B = A/kVolT2M;

value = disc.*(kVolT2M*normpdf(B) - A.*(1-normcdf(B)));
delta = disc.*(1-normcdf(B));
end
