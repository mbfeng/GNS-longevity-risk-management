function [value, delta] = qcallspread_exact(kTau, a, b, x, T, tau, r, AP, EP, kmu, kVol)
% QCALLSPREAD_EXACT  Value and delta at tau of a q-call-spread by numerical integration.
%   Given kappa_tau, kappa_T ~ N(kappa_tau + (T-tau)*kmu, (T-tau)*kVol^2), so
%   V = exp(-r(T-tau)) E[min(max((q - AP)/(EP - AP), 0), 1)] with
%   q = 1 - exp(-exp(a_x + b_x*kappa_T)) is a one-dimensional integral. It is computed
%   with MATLAB's adaptive quadrature (error below 1e-10). Used for the "true" values.
%
% Inputs:
%   kTau  : M-by-1 values of kappa_tau
%   a, b  : Lee-Carter alpha and beta, indexed by age (1 = age 60);  x: reference age index
%   T, tau: maturity and risk horizon, in years;  r: discount rate, continuously compounded
%   AP,EP : attachment and exhaustion points
%   kmu, kVol : drift and volatility of kappa per year
% Outputs:
%   value, delta : M-by-1 value at tau and its derivative in kappa_tau

T2M = T - tau;
disc = exp(-r*T2M);
sd = kVol*sqrt(T2M);
ax = a(x);  bx = b(x);
kAP = (log(-log(1 - AP)) - ax) / bx;      % kappa_T at which q = AP
kEP = (log(-log(1 - EP)) - ax) / bx;      % kappa_T at which q = EP (bx > 0: q increases in kappa)
q  = @(k) 1 - exp(-exp(ax + bx*k));
dq = @(k) bx * exp(ax + bx*k) .* exp(-exp(ax + bx*k));

value = zeros(size(kTau));
delta = zeros(size(kTau));
for i = 1:numel(kTau)
    mu = kTau(i) + T2M*kmu;
    dens = @(k) normpdf(k, mu, sd);
    % payoff is linear in q between the kinks and 1 above EP
    value(i) = disc * (integral(@(k) (q(k) - AP) / (EP - AP) .* dens(k), kAP, kEP, ...
        'AbsTol', 1e-12, 'RelTol', 1e-10) + normcdf(kEP, mu, sd, 'upper'));
    delta(i) = disc * integral(@(k) dq(k) / (EP - AP) .* dens(k), kAP, kEP, ...
        'AbsTol', 1e-12, 'RelTol', 1e-10);
end
end
