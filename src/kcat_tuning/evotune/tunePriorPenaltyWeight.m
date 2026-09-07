function sweepTable = tunePriorPenaltyWeight(ecModel,varargin)
% tunePriorPenaltyWeight  Sweep priorPenaltyWeight and report the
% fit/reproducibility trade-off, so the value doesn't have to be picked
% blind.
%
% priorPenaltyWeight charges for moving a kcat away from its prior, and
% mainly exists to keep large corrections reproducible across seeds
% rather than landing on an arbitrary point along a flat direction. How
% strong a charge is enough is model- and data-specific: too little and
% the search's largest corrections disagree across seeds; too much and
% the search stops finding real corrections at all. This runs
% cmaesKcatTuning at every value in candidates, each at every seed in
% seeds, all against the same tunable set, and reports both fit and
% cross-seed reproducibility per value.
%
% This is exactly as expensive as it sounds: numel(candidates)*numel(seeds)
% full cmaesKcatTuning runs, none of them reusable across candidates. Pass
% a screen you already have to at least avoid recomputing that part.
%
% Parameters
% ----------
% ecModel : struct
%     ecModel with a populated ec.kcat. Not mutated (MATLAB structs are
%     copied on assignment; every candidate/seed pair starts fresh from
%     this input regardless of what earlier runs did).
%
% Name-Value Arguments
% --------------------
% modelAdapter, evoData, tunableMask, screen, targetImpactShare, popsize,
% nProc, verbose
%     as in cmaesKcatTuning.
% candidates : double
%     priorPenaltyWeight values to try (default [0.0, 0.01, 0.03, 0.1]).
% seeds : double
%     seeds to run each candidate at (default [0, 1]). Reproducibility
%     columns pool every pair of seeds, so more than two sharpens them
%     but costs proportionally more.
% foldThreshold : double
%     a kcat counts as "moved" past this fold change from its prior
%     (default 2.0). Matches nChanged's convention loosely (a different
%     quantity, same spirit).
%
% Returns
% -------
% sweepTable : table
%     one row per candidate: priorPenaltyWeight, nSeeds, distanceMean,
%     distanceSpread (max minus min across seeds), nChangedMean,
%     nMovers (kcats past foldThreshold in either seed of a pair, summed
%     over every pair), nBothMoved (past it in both), pctDirectionAgree,
%     medianFoldSpread and maxFoldSpread (among movers, the largest
%     fold-change between what two seeds landed on), and fitCost
%     (distanceMean relative to the best candidate's, as a fraction).
%
% See also
% --------
% cmaesKcatTuning, screenKcatLeverage, selectTunableMask, nChanged
p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'evoData', []; ...
    'tunableMask', []; ...
    'screen', []; ...
    'targetImpactShare', 0.9; ...
    'popsize', []; ...
    'nProc', 1; ...
    'verbose', true; ...
    'candidates', [0.0, 0.01, 0.03, 0.1]; ...
    'seeds', [0, 1]; ...
    'foldThreshold', 2.0});
modelAdapter      = p.modelAdapter;
evoData           = p.evoData;
tunableMask       = p.tunableMask;
screen            = p.screen;
targetImpactShare = p.targetImpactShare;
popsize           = p.popsize;
nProc             = p.nProc;
verbose           = p.verbose;
candidates        = p.candidates;
seeds             = p.seeds;
foldThreshold     = p.foldThreshold;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end
if isempty(evoData)
    evoData = loadEvotuneData('modelAdapter', modelAdapter);
end
if isempty(tunableMask)
    if isempty(screen)
        screen = screenKcatLeverage(ecModel, 'modelAdapter',modelAdapter, ...
            'evoData',evoData, 'nProc',nProc);
    end
    tunableMask = selectTunableMask(ecModel, screen, 'targetImpactShare',targetImpactShare);
end

nCand = numel(candidates);
priorPenaltyWeight = zeros(nCand,1);
nSeeds             = zeros(nCand,1);
distanceMean        = zeros(nCand,1);
distanceSpread       = zeros(nCand,1);
nChangedMean        = zeros(nCand,1);
nMovers              = zeros(nCand,1);
nBothMoved           = zeros(nCand,1);
pctDirectionAgree   = zeros(nCand,1);
medianFoldSpread    = zeros(nCand,1);
maxFoldSpread        = zeros(nCand,1);

for c = 1:nCand
    lam = candidates(c);
    lamParams = modelAdapter.params.evotune;
    lamParams.priorPenaltyWeight = lam;

    vectors = cell(numel(seeds),1);
    distances = zeros(numel(seeds),1);
    changed = zeros(numel(seeds),1);
    oldKcat = [];
    for s = 1:numel(seeds)
        seedVal = seeds(s);
        adapterCopy = modelAdapter;
        adapterCopy.params.evotune = lamParams;
        [~,res] = cmaesKcatTuning(ecModel, 'modelAdapter',adapterCopy, ...
            'evoData',evoData, 'tunableMask',tunableMask, ...
            'popsize',popsize, 'nProc',nProc, 'seed',seedVal, 'verbose',verbose);
        vectors{s} = res.newKcat;
        if ~isempty(res.rmseTrace)
            distances(s) = res.rmseTrace(end);
        else
            distances(s) = NaN;
        end
        changed(s) = nChanged(res.newKcat, res.oldKcat);
        oldKcat = res.oldKcat;
        if verbose
            fprintf('priorPenaltyWeight=%g seed=%d: distance %.4f, %d changed\n', ...
                lam, seedVal, distances(s), changed(s));
        end
    end

    movers = 0; both = 0; agree = 0; spreads = [];
    for i = 1:numel(vectors)
        for j = (i+1):numel(vectors)
            fa = foldChangeLocal(vectors{i}, oldKcat);
            fb = foldChangeLocal(vectors{j}, oldKcat);
            moversIJ = (fa > foldThreshold) | (fb > foldThreshold);
            bothIJ   = (fa > foldThreshold) & (fb > foldThreshold);
            agreeIJ  = sign(log(vectors{i}(bothIJ) ./ oldKcat(bothIJ))) == ...
                       sign(log(vectors{j}(bothIJ) ./ oldKcat(bothIJ)));
            movers = movers + sum(moversIJ);
            both   = both + sum(bothIJ);
            agree  = agree + sum(agreeIJ);
            if any(moversIJ)
                ratio = vectors{i}(moversIJ) ./ vectors{j}(moversIJ);
                spreads = [spreads; max(ratio, 1./ratio)]; %#ok<AGROW>
            end
        end
    end

    priorPenaltyWeight(c) = lam;
    nSeeds(c)              = numel(seeds);
    distanceMean(c)        = mean(distances);
    if numel(distances) > 1
        distanceSpread(c) = max(distances) - min(distances);
    else
        distanceSpread(c) = 0;
    end
    nChangedMean(c) = mean(changed);
    nMovers(c)      = movers;
    nBothMoved(c)   = both;
    if both > 0
        pctDirectionAgree(c) = agree / both;
    else
        pctDirectionAgree(c) = NaN;
    end
    if ~isempty(spreads)
        medianFoldSpread(c) = median(spreads);
        maxFoldSpread(c)    = max(spreads);
    else
        medianFoldSpread(c) = NaN;
        maxFoldSpread(c)    = NaN;
    end
end

fitCost = distanceMean / min(distanceMean) - 1.0;

sweepTable = table(priorPenaltyWeight,nSeeds,distanceMean,distanceSpread,nChangedMean, ...
    nMovers,nBothMoved,pctDirectionAgree,medianFoldSpread,maxFoldSpread,fitCost);
end

function f = foldChangeLocal(kcat,kcat0)
% Per-kcat fold change from the prior, always >= 1. The only caller of
% this in the whole port is this file, so it stays local rather than
% getting its own file (see nChanged.m for the parsimony function that
% genuinely is needed from two places).
f = exp(abs(log(kcat(:) ./ kcat0(:))));
end
