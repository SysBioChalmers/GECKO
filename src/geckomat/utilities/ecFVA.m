function [minFlux, maxFlux] = ecFVA(ecModel, model, varargin)
% ecFVA  Run flux variability analysis on an ecModel.
%
% Flux variability analysis is performed on the ecModel, and isozymic
% reactions are combined to construct output minFlux and maxFlux vectors,
% which follow the same order of model.rxns. The output from this function
% does not include enzyme usage reactions; to observe these, one could
% consider running flux variability directly on the ecModel.
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% model : struct
%     non-ecModel variant of the ecModel, to which the minFlux and maxFlux
%     will be mapped.
%
% Name-Value Arguments
% --------------------
% runParallel : logical
%     speed up calculations by parallel processing. Requires the MATLAB
%     Parallel Computing Toolbox. If this is not installed, the
%     calculations will not be parallelized, regardless of what is
%     indicated as runParallel (default true).
%
% Returns
% -------
% minFlux : double
%     vector of minimum flux rates, corresponding to model.rxns.
% maxFlux : double
%     vector of maximum flux rates, corresponding to model.rxns.
%
% Examples
% --------
%     [minFlux, maxFlux] = ecFVA(ecModel, model);
%
% See also
% --------
% mapRxnsToConv, plotEcFVA

p = parseGECKOargs(varargin, { ...
    'runParallel', true});
runParallel = p.runParallel;

rxnIDs = regexprep(ecModel.rxns,'(_REV)?(_EXP_\d+)?','');
[rxnIDmap, convRxnID] = findgroups(rxnIDs);

N = numel(convRxnID);
maxFluxByGroup = nan(N,1);
minFluxByGroup = nan(N,1);

% nW == 0 forces parfor to run as a plain for loop (no pool required or
% created); nW == Inf lets Parallel Computing Toolbox manage the pool.
% Falls back to serial automatically when the toolbox isn't installed,
% rather than erroring the way an unconditional gcp/parpool call would.
nW = parallelWorkersRAVEN(runParallel);

PB = progressReport(N, 'Running ecFVA');

parfor (i=1:N, nW)
    tmpModel = ecModel;
    tmpModel.c = zeros(numel(tmpModel.rxns),1);

    rxnsToOptim = find(rxnIDmap == i);
    rxnsOptIDs  = ecModel.rxns(rxnsToOptim);
    rxnsToMin = endsWith(rxnsOptIDs,'_REV') | contains(rxnsOptIDs,'_REV_EXP_');
    rxnsToMax = rxnsToOptim(~rxnsToMin);
    rxnsToMin = rxnsToOptim(rxnsToMin);

    tmpModel.c(rxnsToMax) = 1;
    tmpModel.c(rxnsToMin) = -1;
    solMax=solveLP(tmpModel);
    if ~isempty(solMax.x)
        % Exact bound for this canonical reaction, read directly off this
        % group's own solve (the "diagonal"): forward minus reverse
        % variants from the one LP that jointly optimised them, not the
        % best value each variant happened to reach in whichever group's
        % solve maximised it -- an "envelope" that can combine values from
        % different, not necessarily jointly feasible, flux distributions.
        maxFluxByGroup(i) = sum(solMax.x(rxnsToMax)) - sum(solMax.x(rxnsToMin));
    end
    tmpModel.c(rxnsToMax) = -1;
    tmpModel.c(rxnsToMin) = 1;
    solMin=solveLP(tmpModel);
    if ~isempty(solMin.x)
        minFluxByGroup(i) = sum(solMin.x(rxnsToMax)) - sum(solMin.x(rxnsToMin));
    end
    count(PB);
end
PB.done;

% Reorder from canonical-id order to model.rxns order. No min/max swap
% correction is needed here (unlike a row-reduced envelope): both bounds
% come from the same group's LP over the same feasible region with only
% the objective direction flipped, so minFluxByGroup(i) <= maxFluxByGroup(i)
% always holds whenever both solves succeeded.
[mapCheck, origIdx] = ismember(model.rxns, convRxnID);
if ~all(mapCheck)
    error('Not all reactions from model.rxns can be found in the ecModel. Are you sure that ecModel is derived from model?')
end
minFlux = minFluxByGroup(origIdx);
maxFlux = maxFluxByGroup(origIdx);

end
