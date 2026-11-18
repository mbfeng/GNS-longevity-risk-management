function [valTrue, valAnchor] = true_values_mc(kOuter, khat, kmu, kVol, nSteps, nTrue, seed, valueFcn, approxFcn)
% TRUE_VALUES_MC  "True" values V(kappa_tau^(i)) by Monte Carlo with nTrue paths (Section 4).
%   Used when no closed form exists. Every outer scenario i gets its own independent
%   set of nTrue paths (substream i of the random stream), as in the published results.
%   Independent errors average out in the IMSE, eq. (28). (A single shared set of paths
%   would make the errors one common shift, which biases the IMSE of the biased
%   probit-Taylor approximation by about 3%; see README.md, Notes.)
%
% Inputs:
%   kOuter, khat : M-by-1 outer scenarios and the anchor scenario
%   kmu, kVol    : drift and volatility of kappa per year
%   nSteps       : years from tau to maturity (T - tau)
%   nTrue        : number of paths per scenario
%   seed         : random seed
%   valueFcn     : handle, valueFcn(kProj) returns the N-by-1 discounted payoffs
%   approxFcn    : handle, approxFcn(kProj) returns the M-by-1 (probit-)Taylor
%                  approximation from paths started at khat
% Outputs:
%   valTrue   : M-by-1 estimates of V at the outer scenarios
%   valAnchor : M-by-1 (probit-)Taylor approximation from nTrue paths started at khat
%               (substream M+1), i.e. almost without simulation error (Fig. 9)

M = numel(kOuter);
chunk = 2e5;                        % paths per chunk (memory)
valTrue = zeros(M, 1);
parfor i = 1:M
    stream = RandStream('mrg32k3a', 'Seed', seed);
    stream.Substream = i;
    s = 0;
    for j0 = 1:chunk:nTrue
        n = min(chunk, nTrue - j0 + 1);
        s = s + sum(valueFcn(simulate_kappa_paths(kOuter(i), kmu, kVol, nSteps, n, stream)));
    end
    valTrue(i) = s / nTrue;
end
stream = RandStream('mrg32k3a', 'Seed', seed);
stream.Substream = M + 1;
valAnchor = approxFcn(simulate_kappa_paths(khat, kmu, kVol, nSteps, nTrue, stream));
end
