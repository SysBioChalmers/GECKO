function [ecModel,result] = cmaesKcatTuning(ecModel,varargin)
% cmaesKcatTuning  Tune kcats against experimental data with CMA-ES.
%
% Configuration comes entirely from modelAdapter.params.evotune: trust
% tiers (kcatSources/sigma0logSource/sigma0logDefault), maxGrowthWeight,
% priorPenaltyWeight and tieIsozymes control the search, and
% maxGenerations/rmseThreshold are its stopping conditions.
%
% Vendors a separable (diagonal-covariance) CMA-ES -- MATLAB has no
% built-in implementation and Global Optimization Toolbox doesn't ship
% one either. Separable rather than full-covariance because that is what
% geckopy's own port actually runs (cma's CMA_diagonal=True) and what its
% defaults were validated against; it is also lower risk to implement
% correctly. See the "CMA-ES core" local functions below.
%
% Parameters
% ----------
% ecModel : struct
%     ecModel with a populated ec.kcat.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% evoData : struct
%     as returned by loadEvotuneData (default: loaded from modelAdapter).
% tunableMask, screen, targetImpactShare
%     control which kcats are searched over, in order of precedence. Pass
%     tunableMask to fix the set yourself. Otherwise, a mask is built by
%     selectTunableMask at targetImpactShare (default 0.9), from screen if
%     given, or from a fresh screenKcatLeverage call otherwise.
% popsize : double
%     CMA-ES population size. Defaults to 4 + floor(3*log(n)), the
%     standard dimension-scaled default, so it adapts to however many
%     free parameters the mask/tying produced.
% nProc : double
%     number of parallel workers used to score each generation's
%     candidates (default 1, i.e. serial). Requires Parallel Computing
%     Toolbox for > 1.
% seed : double
%     if given, seeds the random draws for a reproducible run.
% verbose : logical
%     whether per-generation progress is printed (default true).
%
% Returns
% -------
% ecModel : struct
%     mutated: tunable rows carry the best kcat vector CMA-ES found, with
%     kcat constraints already applied. MATLAB structs are copied on
%     assignment, so -- unlike geckopy's in-place mutation -- this output
%     must be captured by the caller.
% result : struct
%     rxns/oldKcat/newKcat/groups cover the selected tunable set, tied
%     groups included (their members share one value in newKcat).
%     rmseTrace/objectiveTrace record the best-so-far plain RMSE and
%     optimised objective per generation. nGenerations is how many
%     generations actually ran. converged is true if
%     rmseTrace(end) <= params.rmseThreshold.
%
% Raises
% ------
% error
%     if there are no tunable kcats, or fewer than two free parameters
%     remain after masking and tying -- too little for CMA-ES to search
%     over.
%
% See also
% --------
% screenKcatLeverage, selectTunableMask, tunePriorPenaltyWeight
p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'evoData', []; ...
    'tunableMask', []; ...
    'screen', []; ...
    'targetImpactShare', 0.9; ...
    'popsize', []; ...
    'nProc', 1; ...
    'seed', []; ...
    'verbose', true});
modelAdapter      = p.modelAdapter;
evoData           = p.evoData;
tunableMask       = p.tunableMask;
screen            = p.screen;
targetImpactShare = p.targetImpactShare;
popsize           = p.popsize;
nProc             = max(1, round(p.nProc));
seed              = p.seed;
verbose           = p.verbose;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params.evotune;
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

[tunableIdx,ecRxnIdsTunable,kcat0,groups,sigma0log,tieMap] = ...
    resolveTunableContext(ecModel,params,tunableMask);

reps = find(tieMap == (1:numel(tieMap))');
if numel(reps) < 2
    error('Only %d free parameter(s) after masking and tying; too few for CMA-ES. Widen tunableMask/targetImpactShare.', numel(reps))
end
[~,assign] = ismember(tieMap, reps);  % assign(i): which column of reps tieMap(i) follows

n = numel(reps);
x0  = log(kcat0(reps));
sig = sigma0log(reps);
[loFull,hiFull] = kcatBoundsLocal(kcat0);
loLog = log(loFull(reps));
hiLog = log(hiFull(reps));

if isempty(popsize)
    popsize = 4 + floor(3*log(n));
end
popsize = max(4, round(popsize));

expand = @(Xlog) exp(Xlog(assign,:));

% Baseline: score the untouched prior, and pick up the excarbon-augmented
% model to carry through the rest of the run (see evotuneScore's 4th
% output -- avoids rescanning metabolite formulas on every candidate).
[obj0,rmse0,~,baseModel] = scoreCandidate(ecModel,tunableIdx,ecRxnIdsTunable,evoData,modelAdapter,params,kcat0,kcat0,sigma0log);
bestObj = obj0; bestRmse = rmse0; bestVec = kcat0;

rmseTrace = []; objectiveTrace = [];

[meanX,sigma,D,pSigma,pC,cs,ds,cc,c1,cmu,muEff,w,mu,chiN] = cmaesInit(x0,sig,popsize);

if ~isempty(seed)
    rng(seed);
end

generation = 0;
while generation < params.maxGenerations
    Z = randn(n,popsize);
    Xlog = meanX + sigma * (D .* Z);
    Xlog = min(max(Xlog,loLog),hiLog);
    % Recompute Z consistent with any clipping (D>0 throughout by construction).
    Z = (Xlog - meanX) ./ (sigma * D);

    kcatMatrix = expand(Xlog);
    objVals = zeros(popsize,1); rmseVals = zeros(popsize,1);
    if nProc == 1
        for k = 1:popsize
            [objVals(k),rmseVals(k),~,baseModel] = scoreCandidate( ...
                baseModel,tunableIdx,ecRxnIdsTunable,evoData,modelAdapter,params, ...
                kcatMatrix(:,k),kcat0,sigma0log);
        end
    else
        parfor k = 1:popsize
            [objVals(k),rmseVals(k)] = scoreCandidate( ...
                ecModel,tunableIdx,ecRxnIdsTunable,evoData,modelAdapter,params, ...
                kcatMatrix(:,k),kcat0,sigma0log);
        end
    end
    generation = generation + 1;

    [~,ord] = sort(objVals,'ascend');
    if objVals(ord(1)) < bestObj
        bestObj = objVals(ord(1));
        bestRmse = rmseVals(ord(1));
        bestVec = kcatMatrix(:,ord(1));
    end
    rmseTrace(end+1) = bestRmse; %#ok<AGROW>
    objectiveTrace(end+1) = bestObj; %#ok<AGROW>
    if verbose
        fprintf('generation %d: objective %.4f, rmse %.4f\n', generation, bestObj, bestRmse);
    end
    if bestRmse <= params.rmseThreshold
        break
    end

    [meanX,sigma,D,pSigma,pC] = cmaesUpdate( ...
        meanX,sigma,D,pSigma,pC,Z,ord,w,mu,cs,ds,cc,c1,cmu,muEff,chiN,generation);
end

ecModel.ec.kcat(tunableIdx) = bestVec;
ecModel = applyKcatConstraints(ecModel,'updateRxns',ecRxnIdsTunable);

result = struct();
result.rxns            = ecRxnIdsTunable;
result.oldKcat          = kcat0;
result.newKcat          = bestVec;
result.groups           = groups;
result.rmseTrace        = rmseTrace;
result.objectiveTrace   = objectiveTrace;
result.nGenerations     = generation;
result.converged        = bestRmse <= params.rmseThreshold;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scoring core, shared with screenKcatLeverage/tunePriorPenaltyWeight.

function [obj,rmse,rmseDetail,ecModel] = scoreCandidate(ecModel,tunableIdx,ecRxnIdsTunable,evoData,modelAdapter,params,kcatVec,kcat0,sigma0log)
ecModel.ec.kcat(tunableIdx) = kcatVec;
ecModel = applyKcatConstraints(ecModel,'updateRxns',ecRxnIdsTunable);
[obj,rmse,rmseDetail,ecModel] = evotuneScore(ecModel,evoData, ...
    'modelAdapter',modelAdapter, 'maxGrowthWeight',params.maxGrowthWeight, ...
    'priorPenaltyWeight',params.priorPenaltyWeight, ...
    'kcatVec',kcatVec, 'kcat0',kcat0, 'sigma0log',sigma0log);
end

function [lo,hi] = kcatBoundsLocal(kcat0)
% Biologically plausible bounds for proposed kcats, in 1/s. The window is
% 1e-2 to 1e4 for an ordinary kcat, and always widened far enough to
% contain the prior with a hundred-fold margin either side.
lo = min(1e-2, kcat0/100);
hi = max(1e4, kcat0*100);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Separable (diagonal) CMA-ES core.
% Standard (mu/mu_w,lambda)-CMA-ES with a diagonal-only covariance update
% (Ros & Hansen 2008), operating in log-kcat space. Sampling, ranking,
% recombination, step-size (cumulative step-size adaptation) and
% per-coordinate variance updates are the textbook equations; see e.g.
% Hansen, "The CMA Evolution Strategy: A Tutorial".

function [meanX,sigma,D,pSigma,pC,cs,ds,cc,c1,cmu,muEff,w,mu,chiN] = cmaesInit(x0,sig,popsize)
n = numel(x0);
mu = floor(popsize/2);
wRaw = log(mu+0.5) - log((1:mu)');
w = wRaw / sum(wRaw);
muEff = 1 / sum(w.^2);

cs  = (muEff+2) / (n+muEff+5);
ds  = 1 + 2*max(0,sqrt((muEff-1)/(n+1))-1) + cs;
cc  = (4+muEff/n) / (n+4+2*muEff/n);
c1  = 2 / ((n+1.3)^2+muEff);
cmu = min(1-c1, 2*(muEff-2+1/muEff) / ((n+2)^2+muEff));
chiN = sqrt(n) * (1 - 1/(4*n) + 1/(21*n^2));

meanX = x0;
% Initial per-coordinate step scaled by each parameter's own prior sigma
% (relative to the median), rather than uniform -- mirrors passing
% scaling_of_variables to cma in the Python port.
D = sig / median(sig);
sigma = median(sig);
pSigma = zeros(n,1);
pC = zeros(n,1);
end

function [meanX,sigma,D,pSigma,pC] = cmaesUpdate( ...
    meanX,sigma,D,pSigma,pC,Z,ord,w,mu,cs,ds,cc,c1,cmu,muEff,chiN,generation)
n = numel(meanX);
Zbest = Z(:,ord(1:mu));
zw = Zbest * w;

meanX = meanX + sigma * (D .* zw);

pSigma = (1-cs)*pSigma + sqrt(cs*(2-cs)*muEff) * zw;
sigma = sigma * exp((cs/ds) * (norm(pSigma)/chiN - 1));

hSigmaThresh = (1.4 + 2/(n+1)) * chiN;
hSigma = (norm(pSigma) / sqrt(1-(1-cs)^(2*generation))) < hSigmaThresh;

pC = (1-cc)*pC + hSigma * sqrt(cc*(2-cc)*muEff) * (D .* zw);

Cdiag = D.^2;
rankMuTerm = (D.*Zbest).^2 * w;
Cdiag = (1-c1-cmu)*Cdiag + c1*(pC.^2 + (1-hSigma)*cc*(2-cc)*Cdiag) + cmu*rankMuTerm;
Cdiag = max(Cdiag, 1e-300);
D = sqrt(Cdiag);
end
