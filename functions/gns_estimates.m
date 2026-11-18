function [vGns, vGnsSn] = gns_estimates(kOuter, kmu, kVol, val, samp)
% GNS_ESTIMATES  Green nested simulation (GNS) estimates of V(kappa_tau) for all outer scenarios.
%   The N inner paths were simulated in equal numbers from every outer scenario, so
%   their sampling density is the M-component mixture of eq. (26),
%   g(s) = mean_i f(s | kOuter_i), where f is the one-year normal density of kappa.
%   Every path is reused for every scenario with the likelihood ratio
%   LR_j(i) = f(samp_j | kOuter_i) / g(samp_j) (Section 3.1).
%
% Inputs:
%   kOuter : M-by-1 outer scenarios kappa_tau^(i)
%   kmu    : drift of kappa per year
%   kVol   : volatility of kappa per year
%   val    : N-by-1 discounted payoffs of the inner paths
%   samp   : N-by-1 kappa one year after tau on each inner path
% Outputs:
%   vGns   : M-by-1 regular GNS estimates, mean(val .* LR)
%   vGnsSn : M-by-1 self-normalized GNS estimates, mean(val .* LR) / mean(LR)

M = numel(kOuter);
N = numel(samp);
c = 1 / (kVol*sqrt(2*pi));
s1 = zeros(M, 1);                     % sum over paths of LR .* val
s0 = zeros(M, 1);                     % sum over paths of LR
block = max(1, floor(1e7 / M));       % paths per block: each M-by-n matrix is about 80 MB
                                      % (a few such matrices exist at once, per worker)
for j0 = 1:block:N
    jj = j0:min(j0 + block - 1, N);
    F = c * exp(-0.5 * ((samp(jj)' - kOuter(:) - kmu) / kVol).^2);   % f(samp_j | kOuter_i), M-by-n
    LR = F ./ mean(F, 1);                                             % divide by mixture density g
    s1 = s1 + LR * val(jj);
    s0 = s0 + sum(LR, 2);
end
vGns = s1 / N;
vGnsSn = s1 ./ s0;
end
