function lc = kappa_model()
% KAPPA_MODEL  Lee-Carter fit and the random walk with drift for kappa (Section 4).
%   Fits the Lee-Carter model to data/hmd_ew_female_1962_2011.mat (England & Wales
%   females, ages 60-89, 1962-2011) by Poisson maximum likelihood, then fits
%   diff(k) ~ N(kmu, kVar) by maximum likelihood (sample mean and variance with 1/n).
%
% Output: struct with fields
%   a, b  : 30-by-1 alpha_x and beta_x (index 1 = age 60)
%   k     : 50-by-1 kappa_t (1962-2011)
%   k0    : kappa in 2011 (time 0)
%   kmu, kVar, kVol : drift, variance and volatility of kappa per year
%   Age, Year : the ages and years of the data

p = project_paths();
H = load(fullfile(p.data, 'hmd_ew_female_1962_2011.mat'));   % Age, Year, Death, Expos
[a, b, k] = lc_fit_svd(H.Expos, H.Death);
[a, b, k] = lc_fit_mle(H.Expos, H.Death, a, b, k);

dk = diff(k);
lc = struct('a', a, 'b', b, 'k', k, 'k0', k(end), 'kmu', mean(dk), 'kVar', var(dk, 1), ...
    'kVol', sqrt(var(dk, 1)), 'Age', H.Age, 'Year', H.Year);
end
