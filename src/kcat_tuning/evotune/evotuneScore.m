function [objective,rmse,rmseDetail,ecModel] = evotuneScore(ecModel,evoData,varargin)
% evotuneScore  Score one ecModel's current kcats against evotune data.
%
% Internal helper used by cmaesKcatTuning, screenKcatLeverage and
% tunePriorPenaltyWeight to score one ecModel (with a given set of kcats
% already applied via applyKcatConstraints) against the experimental
% growth and flux data. Fuses what geckopy splits into simulate.py (the
% FBA half) and distance.py (the RMSE half) into one function, matching
% this module's own predecessor abc_max.m -- there is no MATLAB test
% suite here to benefit from Python's testability-driven split.
%
% Parameters
% ----------
% ecModel : struct
%     ecModel with the current candidate kcats already applied (via
%     applyKcatConstraints). Not mutated.
% evoData : struct
%     structure with experimental data to be used, as returned by
%     loadEvotuneData.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% maxGrowthWeight : double
%     relative weight of the max-growth dataset against the flux dataset:
%     (rmseFlux + w*rmseMaxGrate) / (w+1) (default 1.0).
% priorPenaltyWeight : double
%     weight on a prior term added to the objective (default 0, i.e. the
%     objective equals the plain RMSE). See kcatVec/kcat0/sigma0log below.
% kcatVec, kcat0, sigma0log : double
%     only used when priorPenaltyWeight > 0: the tunable kcat subset's
%     current values, prior values, and per-kcat prior log-space std dev,
%     all the same length and in the same order. The penalty added to the
%     objective is priorPenaltyWeight * mean((log(kcatVec./kcat0)./sigma0log).^2).
%
% Returns
% -------
% objective : double
%     rmse, plus the prior penalty term when priorPenaltyWeight > 0. This
%     is what a search should minimise.
% rmse : double
%     the plain (unpenalised) average RMSE, always comparable across runs
%     regardless of priorPenaltyWeight.
% rmseDetail : cell
%     per-condition RMSE values, labelled by data type and condition index.
% ecModel : struct
%     the input ecModel, with ec.excarbon cached on it if it wasn't
%     already present. MATLAB structs are copied on assignment, not
%     shared by reference like a Python object, so callers scoring many
%     candidates in a loop should carry this output back in as the base
%     copy for the next candidate, rather than recomputing excarbon (an
%     expensive metabolite-formula scan) on every single call.
%
% See also
% --------
% cmaesKcatTuning, screenKcatLeverage, loadEvotuneData

p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'maxGrowthWeight', 1.0; ...
    'priorPenaltyWeight', 0.0; ...
    'kcatVec', []; ...
    'kcat0', []; ...
    'sigma0log', []});
modelAdapter        = p.modelAdapter;
maxGrowthWeight     = p.maxGrowthWeight;
priorPenaltyWeight  = p.priorPenaltyWeight;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

if ~isfield(ecModel,'excarbon')
    ecModel = fillCarbonNum(ecModel);
    ecModel.excarbon(ecModel.excarbon == 0) = 1;
end

rmse_1 = []; rmseDetail1 = []; rmse_2 = []; rmseDetail2 = [];

if ~isempty(evoData.fluxData)
    [rmse_1, rmseDetail1] = simulateCondition(ecModel,evoData.fluxData,true,evoData.zeroFlux,modelAdapter);
    rmseDetail1 = [cellstr("fluxData_" + string(1:numel(rmseDetail1)))',num2cell(rmseDetail1)];
end
if ~isempty(evoData.maxGrate)
    [rmse_2, rmseDetail2] = simulateCondition(ecModel,evoData.maxGrate,false,evoData.zeroFlux,modelAdapter);
    rmseDetail2 = [cellstr("maxGrowth_" + string(1:numel(rmseDetail2)))',num2cell(rmseDetail2)];
end

parts = []; weights = [];
if ~isempty(rmse_1)
    parts(end+1) = rmse_1; weights(end+1) = 1.0;
end
if ~isempty(rmse_2)
    parts(end+1) = rmse_2; weights(end+1) = maxGrowthWeight;
end
if isempty(parts)
    rmse = NaN;
else
    rmse = sum(parts .* weights) / sum(weights);
end
rmseDetail = [rmseDetail1;rmseDetail2];

objective = rmse;
if priorPenaltyWeight > 0 && ~isempty(p.kcatVec) && ~isempty(p.kcat0) && ~isempty(p.sigma0log)
    dev = log(p.kcatVec(:) ./ p.kcat0(:)) ./ p.sigma0log(:);
    objective = rmse + priorPenaltyWeight * mean(dev.^2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rmse, rmseList] = simulateCondition(ecModel,data,constrain,rxn2block,modelAdapter)
% One condition's FBA simulation + carbon-weighted RMSE against measured
% data. Ported from abc_max.m's rmsecal.

rmseList = zeros(length(data.conds),1);

[~,exchIdx] = ismember(data.conds,data.exchMets);
if any(exchIdx == 0)
    missingMets = strjoin(unique(data.conds(exchIdx == 0)), '; ');
    error('Carbon source(s) "%s" in the provided fluxData or maxGrowth cannot be matched by name with an exchange reaction.', missingMets);
end

for i = 1:length(data.conds)
    % Set all other carbon sources to zero, then unblock this condition's
    model_tmp = setParam(ecModel,'lb',data.exchRxnIDs(unique(exchIdx)),0);
    if constrain % RMSE from flux data: constrain carbon uptake at the measured rate
        model_tmp = setParam(model_tmp,'lb',data.exchRxnIDs(exchIdx(i)),data.exchFluxes(i,exchIdx(i)));
    else % RMSE from max growth: leave carbon uptake fully open
        model_tmp = setParam(model_tmp,'lb',data.exchRxnIDs(exchIdx(i)),-1000);
    end

    % This row's own carbon source is never itself asserted zero-flux.
    rxn2blockIter = rxn2block;
    rxn2blockIter(strcmp(rxn2blockIter,data.exchRxnIDs(exchIdx(i)))) = [];

    o2Flux = data.exchFluxes(i,strcmp('oxygen',data.exchMets));
    if o2Flux == 0
        model_tmp = modelAdapter.makeModelAnaerobic(model_tmp);
    end
    if ~(isnan(data.Ptot(i)) || data.Ptot(i) == 0)
        model_tmp = modelAdapter.changeProteinBiomass(model_tmp,data.Ptot(i));
    end

    sol = solveLP(model_tmp);
    if checkSolution(sol)
        bioRxn  = getIndexes(model_tmp,data.biomass,'rxns');
        % Biomass carries an assumed 41 Cmmol/gDCW, matching fillCarbonNum.
        bioMeas = data.grRate(i) * 41;
        bioSim  = sol.x(bioRxn) * 41;

        if constrain
            fluxToCheck   = ~isnan(data.exchFluxes(i,:));
            fluxRxnIdx    = getIndexes(model_tmp,data.exchRxnIDs(fluxToCheck),'rxns');
            cNormMeasFlux = model_tmp.excarbon(fluxRxnIdx) .* transpose(data.exchFluxes(i,fluxToCheck));
            cNormSimFlux  = model_tmp.excarbon(fluxRxnIdx) .* sol.x(fluxRxnIdx);

            blockRxnIdx = getIndexes(model_tmp,rxn2blockIter,'rxns');
            blockSim    = sol.x(blockRxnIdx) .* model_tmp.excarbon(blockRxnIdx);
            % Already-correct zero predictions don't dilute the RMSE.
            blockSim(blockSim == 0) = [];
            blockMeas = zeros(numel(blockSim),1);

            measured  = [cNormMeasFlux; bioMeas; blockMeas];
            simulated = [cNormSimFlux; bioSim; blockSim];
        else
            measured  = bioMeas;
            simulated = bioSim;
        end
        rmseList(i) = sqrt(mean((measured - simulated).^2));
    else
        rmseList(i) = NaN;
    end
end
rmseList(isnan(rmseList)) = 99; % Fixed penalty for an infeasible condition.
if isfield(data,'evotuneRMSEweight') && ~isempty(data.evotuneRMSEweight)
    rmseList = rmseList .* data.evotuneRMSEweight;
end
rmse = mean(rmseList);
end
