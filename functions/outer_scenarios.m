function [kOuter, khat, outerPdf] = outer_scenarios(lc, tau, M)
% OUTER_SCENARIOS  Outer scenarios kappa_tau^(i), i = 1..M (Section 4).
%   Instead of random draws, the scenarios are M equally spaced quantiles of
%   kappa_tau ~ N(k0 + tau*kmu, tau*kVar), for better uniformity.
%
% Inputs:  lc (from kappa_model), risk horizon tau in years, number of scenarios M
% Outputs: kOuter   : M-by-1 sorted outer scenarios
%          khat     : anchor scenario of the Taylor approximations (mean of kOuter)
%          outerPdf : M-by-1 density of kappa_tau at kOuter (for the setting figures)

mu = lc.k0 + tau*lc.kmu;
vol = sqrt(tau*lc.kVar);

pct = linspace(0, 1, M + 2)';
pct([1, end]) = [];
kOuter = sort(norminv(pct, mu, vol));
khat = mean(kOuter);
outerPdf = normpdf(kOuter, mu, vol);
end
