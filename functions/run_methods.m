function Val = run_methods(kOuter, khat, kmu, kVol, nSteps, budgets, nMac, seed, valueFcn, approxFcn)
% RUN_METHODS  Repeated estimates of V(kappa_tau) by the four valuation methods (Section 4).
%   For each of nMac macro replications and each simulation budget:
%     Nested : standard nested simulation, budget/M inner paths per outer scenario
%     GNS    : green nested simulation; reuses the same inner paths for every outer
%              scenario with the mixture likelihood ratio, eq. (26)
%     GNS_SN : self-normalized GNS (same paths)
%     Approx : (probit-)Taylor approximation from budget paths started at khat
%   Replications run in parallel (parfor) if the Parallel Computing Toolbox is
%   available. Replication m always uses substream m of the random stream, so the
%   random numbers do not depend on the number of workers (results agree up to
%   floating-point rounding, about 1e-13).
%
% Inputs:
%   kOuter, khat : M-by-1 outer scenarios and the anchor scenario
%   kmu, kVol    : drift and volatility of kappa per year
%   nSteps       : years from tau to maturity (T - tau)
%   budgets      : 1-by-B simulation budgets (total inner paths), e.g. [1e3 1e4 1e5]
%   nMac         : number of macro replications R, eq. (28)
%   seed         : random seed
%   valueFcn     : handle, valueFcn(kProj) returns N-by-1 discounted payoffs
%   approxFcn    : handle, approxFcn(kProj) returns the M-by-1 approximation built
%                  from paths kProj started at khat
% Output:
%   Val : struct with fields Nested, Approx, GNS, GNS_SN, each M-by-nMac-by-B

M = numel(kOuter);
B = numel(budgets);
nIn = budgets / M;                  % inner paths per outer scenario
[nested, approx, gns, gnsSn] = deal(zeros(M, nMac, B));

t0 = tic;
parfor m = 1:nMac
    stream = RandStream('mrg32k3a', 'Seed', seed);
    stream.Substream = m;
    [vN, vA, vG, vS] = deal(zeros(M, B));
    for b = 1:B
        % inner paths: nIn(b) from each outer scenario, shared by nested simulation and GNS
        kProj = simulate_kappa_paths(kOuter, kmu, kVol, nSteps, budgets(b), stream);
        val = valueFcn(kProj);
        vN(:, b) = mean(reshape(val, nIn(b), M), 1)';
        [vG(:, b), vS(:, b)] = gns_estimates(kOuter, kmu, kVol, val, kProj(1, :)');
        % Taylor approximation from paths started at the anchor scenario
        vA(:, b) = approxFcn(simulate_kappa_paths(khat, kmu, kVol, nSteps, budgets(b), stream));
    end
    nested(:, m, :) = reshape(vN, M, 1, B);
    approx(:, m, :) = reshape(vA, M, 1, B);
    gns(:, m, :) = reshape(vG, M, 1, B);
    gnsSn(:, m, :) = reshape(vS, M, 1, B);
end
fprintf('%d macro replications done in %.1f minutes\n', nMac, toc(t0)/60);
Val = struct('Nested', nested, 'Approx', approx, 'GNS', gns, 'GNS_SN', gnsSn);
end
