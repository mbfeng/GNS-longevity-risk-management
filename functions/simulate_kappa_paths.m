function kProj = simulate_kappa_paths(kStart, kmu, kVol, nSteps, nPaths, stream)
% SIMULATE_KAPPA_PATHS  Simulate paths of the Lee-Carter period index kappa.
%   kappa_{t+1} = kappa_t + kmu + kVol*Z, Z ~ N(0,1) (random walk with drift).
%   nPaths/numel(kStart) paths start from each value in kStart; the columns are
%   grouped by starting value (first all paths from kStart(1), then kStart(2), ...).
%
% Inputs:
%   kStart : scalar or vector of starting values (e.g. the outer scenarios kappa_tau)
%   kmu    : drift of kappa per year
%   kVol   : volatility of kappa per year
%   nSteps : number of years to simulate (T - tau for inner paths)
%   nPaths : total number of paths; a multiple of numel(kStart)
%   stream : RandStream to draw from
% Output:
%   kProj  : nSteps-by-nPaths matrix; row t is kappa at tau + t

nReps = nPaths / numel(kStart);
Z = randn(stream, nSteps, nPaths);
kProj = cumsum(kmu + kVol*Z, 1) + repelem(kStart(:)', 1, nReps);
end
