function screen = screenKcatLeverage(ecModel,varargin)
% screenKcatLeverage  Report which kcats the data can actually speak to.
%
% No optimisation happens here -- this only says which parameters are
% worth curating or tuning, and can be read on its own well before any
% tuning run. Each tie-group (isozyme copies sharing a prior and a source
% move together when params.tieIsozymes) is perturbed up and down by
% fold, and its leverage is the largest resulting change in RMSE. Groups
% are ranked by that leverage weighted by sigma0log: a trusted source
% needs proportionally more measured effect to rank alongside an
% untrusted one.
%
% Costs one simulation per condition per tie-group probed (two probes
% each, up and down) plus one baseline -- comparable in scale to a full
% tuning run's evaluation budget, not a cheap pre-check.
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
% tunableMask : logical
%     restrict which kcats are screened (default: every ec.kcat > 0).
% fold : double
%     fold-change used to perturb each tie-group up and down (default 2.0).
% nProc : double
%     number of parallel workers used to score candidate vectors (default
%     1, i.e. serial). Requires Parallel Computing Toolbox for > 1.
%
% Returns
% -------
% screen : table
%     one row per tie-group, sorted by rank descending: rxnId (the
%     group's representative), nIsozymes, sourceGroup, kcat0, leverage,
%     sigma0log, rank (leverage*sigma0log), cumLeverageShare (running
%     total of rank as a fraction of its sum), and positions (a cell
%     array of indices into ecModel.ec.rxns) -- what selectTunableMask
%     needs to expand a row back into a mask.
%
% See also
% --------
% selectTunableMask, cmaesKcatTuning
p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'evoData', []; ...
    'tunableMask', []; ...
    'fold', 2.0; ...
    'nProc', 1});
modelAdapter = p.modelAdapter;
evoData      = p.evoData;
tunableMask  = p.tunableMask;
fold         = p.fold;
nProc        = max(1, round(p.nProc));

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

[tunableIdx,ecRxnIdsTunable,kcat0,groups,sigma0log,tieMap] = ...
    resolveTunableContext(ecModel,params,tunableMask);

reps = find(tieMap == (1:numel(tieMap))');
nReps = numel(reps);
membersOf = cell(nReps,1);
for r = 1:nReps
    membersOf{r} = find(tieMap == reps(r));
end

vectors = cell(1 + 2*nReps, 1);
vectors{1} = kcat0;
for r = 1:nReps
    members = membersOf{r};
    up = kcat0; up(members) = up(members) * fold;
    dn = kcat0; dn(members) = dn(members) / fold;
    vectors{1 + 2*r - 1} = up;
    vectors{1 + 2*r}     = dn;
end
nVec = numel(vectors);
rmses = zeros(nVec,1);

if nProc == 1
    % Serial path: carry the excarbon-augmented model from one call into
    % the next, so the metabolite-formula scan behind it (evotuneScore's
    % 4th output) runs once for the whole screen rather than once per
    % probe vector.
    baseModel = ecModel;
    for v = 1:nVec
        [~,rmses(v),~,baseModel] = scoreVector(baseModel,tunableIdx,ecRxnIdsTunable,evoData,modelAdapter,params,vectors{v});
    end
else
    parfor v = 1:nVec
        [~,rmses(v)] = scoreVector(ecModel,tunableIdx,ecRxnIdsTunable,evoData,modelAdapter,params,vectors{v});
    end
end

baseRmse = rmses(1);
rxnId       = cell(nReps,1);
nIsozymes   = zeros(nReps,1);
sourceGroup = cell(nReps,1);
kcat0Col    = zeros(nReps,1);
leverage    = zeros(nReps,1);
sigma0logCol= zeros(nReps,1);
positions   = cell(nReps,1);
for r = 1:nReps
    upRmse = rmses(1 + 2*r - 1);
    dnRmse = rmses(1 + 2*r);
    members = membersOf{r};
    rxnId{r}        = ecRxnIdsTunable{reps(r)};
    nIsozymes(r)     = numel(members);
    sourceGroup{r}   = groups{reps(r)};
    kcat0Col(r)      = kcat0(reps(r));
    leverage(r)      = max(abs(upRmse - baseRmse), abs(dnRmse - baseRmse));
    sigma0logCol(r)  = sigma0log(reps(r));
    positions{r}     = tunableIdx(members);
end

rnk = leverage .* sigma0logCol;
[~,ord] = sort(rnk,'descend');
rxnId=rxnId(ord); nIsozymes=nIsozymes(ord); sourceGroup=sourceGroup(ord);
kcat0Col=kcat0Col(ord); leverage=leverage(ord); sigma0logCol=sigma0logCol(ord);
positions=positions(ord); rnk=rnk(ord);

total = sum(rnk);
if total > 0
    cumLeverageShare = cumsum(rnk) / total;
else
    cumLeverageShare = zeros(size(rnk));
end

screen = table(rxnId,nIsozymes,sourceGroup,kcat0Col,leverage,sigma0logCol,rnk,cumLeverageShare,positions, ...
    'VariableNames',{'rxnId','nIsozymes','sourceGroup','kcat0','leverage','sigma0log','rank','cumLeverageShare','positions'});
end

function [obj,rmse,rmseDetail,ecModel] = scoreVector(ecModel,tunableIdx,ecRxnIdsTunable,evoData,modelAdapter,params,kcatVec)
ecModel.ec.kcat(tunableIdx) = kcatVec;
ecModel = applyKcatConstraints(ecModel,'updateRxns',ecRxnIdsTunable);
[obj,rmse,rmseDetail,ecModel] = evotuneScore(ecModel,evoData,'modelAdapter',modelAdapter,'maxGrowthWeight',params.maxGrowthWeight);
end
