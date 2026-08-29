function result = relaxProteomicsGreedy(model, varargin)
% relaxProteomicsGreedy  Greedy shadow-price-ordered relaxation of proteomics
% constraints.
%
% While growth is below minimalGrowth, finds the still-proteomics-constrained
% enzyme (usage_prot_<id> upper bound below defaultUpperBound) with the
% largest absolute shadow price on its prot_<id> mass-balance constraint,
% fully relaxes it back to defaultUpperBound, re-solves, and repeats.
%
% Parameters
% ----------
% model : struct
%     a full ecModel in GECKO 3 format, with some usage_prot_<id> reactions
%     constrained below defaultUpperBound (e.g. by constrainEnzConcs).
%
% Name-Value Arguments
% --------------------
% minimalGrowth : double
%     target objective value. Required.
% enzymeSet : cell
%     subset of uniprot ids eligible for relaxation (default: every
%     currently proteomics-constrained enzyme).
% maxIterations : double
%     safety cap; errors if exceeded while eligible enzymes still remained
%     (default 100).
% defaultUpperBound : double
%     value to restore on a relaxed usage reaction, matching the shared-pool
%     default makeEcModel itself sets (default 1000).
%
% Returns
% -------
% result : struct
%     - relaxed : struct mapping each relaxed uniprot id to its original
%       model.ec.concs value, so the caller can restore it later.
%     - trace : struct array, one row per relaxation step, with fields
%       iteration (0-based), relaxedUniprot, growthBefore, growthAfter,
%       shadowPrice.
%     - finalGrowth : objective value when the loop stopped.
%     - converged : true iff finalGrowth >= minimalGrowth.
%
% Examples
% --------
%     result = relaxProteomicsGreedy(ecModel, 'minimalGrowth', 0.1);
%
% Notes
% -----
% Different relaxation strategy from flexibilizeEnzConcs: this one is
% shadow-price-ordered (one solve per step) and often converges in fewer
% iterations when one or two enzymes dominate the infeasibility. The
% trade-off is that it fully unconstrains each picked enzyme with no
% tighten-back pass, so the result is looser (less proteomics-faithful).
%
% There are two distinct non-convergence outcomes: running out of
% eligible candidates before reaching minimalGrowth returns normally with
% converged=false; exhausting maxIterations while candidates still
% remained raises an error instead.
%
% See also
% --------
% flexibilizeEnzConcs, getEnzymeBottlenecks

p = parseGECKOargs(varargin, {'minimalGrowth', []; 'enzymeSet', []; ...
    'maxIterations', []; 'defaultUpperBound', []});
minimalGrowth     = p.minimalGrowth;
enzymeSet         = p.enzymeSet;
maxIterations     = p.maxIterations;
defaultUpperBound = p.defaultUpperBound;
if isempty(minimalGrowth)
    error('relaxProteomicsGreedy: minimalGrowth is required.')
end
if isempty(maxIterations)
    maxIterations = 100;
end
if isempty(defaultUpperBound)
    defaultUpperBound = 1000;
end

if isfield(model.ec, 'geckoLight') && model.ec.geckoLight
    error(['relaxProteomicsGreedy requires a full ecModel, not gecko-light: light models have no ' ...
        'per-enzyme prot_<id> / usage_prot_<id> machinery.'])
end

usageRxns = strcat('usage_prot_', model.ec.enzymes);
[found, usageIdx] = ismember(usageRxns, model.rxns);
if ~all(found)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end
protMets = strcat('prot_', model.ec.enzymes);
[found, protMetIdx] = ismember(protMets, model.mets);
if ~all(found)
    error('Cannot find prot_<id> metabolites for all enzymes. This is done by makeEcModel.')
end

if isempty(enzymeSet)
    pool = model.ec.enzymes;
else
    pool = intersect(enzymeSet(:), model.ec.enzymes, 'stable');
end
[~, poolIdx] = ismember(pool, model.ec.enzymes);
poolIdx = poolIdx(poolIdx > 0);
candidates = poolIdx(model.ub(usageIdx(poolIdx)) < defaultUpperBound);

relaxed = struct();
trace = struct('iteration', {}, 'relaxedUniprot', {}, 'growthBefore', {}, ...
    'growthAfter', {}, 'shadowPrice', {});

sol = solveLP(model);
if sol.stat == 1
    growth = sol.f;
else
    growth = -inf;
end

for it = 0:(maxIterations - 1)
    if growth >= minimalGrowth
        result.relaxed = relaxed;
        result.trace = trace;
        result.finalGrowth = growth;
        result.converged = true;
        return
    end
    if isempty(candidates)
        warning('relaxProteomicsGreedy: no more candidates to relax; stopping.')
        result.relaxed = relaxed;
        result.trace = trace;
        result.finalGrowth = growth;
        result.converged = false;
        return
    end

    shadowPrices = abs(sol.sPrice(protMetIdx(candidates)));
    [~, bestPos] = max(shadowPrices);
    targetIdx = candidates(bestPos);
    targetUniprot = model.ec.enzymes{targetIdx};
    sp = sol.sPrice(protMetIdx(targetIdx));

    relaxed.(matlab.lang.makeValidName(targetUniprot)) = model.ec.concs(targetIdx);
    model.ub(usageIdx(targetIdx)) = defaultUpperBound;
    candidates(candidates == targetIdx) = [];

    newSol = solveLP(model);
    if newSol.stat == 1
        newGrowth = newSol.f;
    else
        newGrowth = -inf;
    end

    trace(end+1) = struct('iteration', it, 'relaxedUniprot', targetUniprot, ...
        'growthBefore', growth, 'growthAfter', newGrowth, 'shadowPrice', sp); %#ok<AGROW>

    growth = newGrowth;
    sol = newSol;
end

error('relaxProteomicsGreedy: did not converge in %d iterations (final growth %.4g < target %.4g).', ...
    maxIterations, growth, minimalGrowth)
end
