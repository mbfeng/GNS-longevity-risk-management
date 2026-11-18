function approx = linear_taylor_sim(kProj, kOuter, khat, valueFcn)
% LINEAR_TAYLOR_SIM  Linear Taylor approximation with simulated value and delta (Section 2.4.3).
%   V(kappa) ~ mean(value) + mean(delta)*(kappa - khat), with value and pathwise delta
%   from paths kProj started at the anchor scenario khat.
%
% Inputs:
%   kProj    : paths started at khat
%   kOuter   : M-by-1 points at which the approximation is evaluated
%   khat     : anchor scenario
%   valueFcn : handle, [value, delta] = valueFcn(kProj)
% Output:
%   approx   : M-by-1 approximate values

[value, delta] = valueFcn(kProj);
approx = mean(value) + mean(delta) * (kOuter - khat);
end
