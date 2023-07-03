function [minFlux, maxFlux] = ecFVA(ecModel, model)
% ecFVA
%   Flux variability analysis is performed on the ecModel, and isozymic
%   reactions are combined to construct ouput minFlux and maxFlux vectors,
%   which follow the same order of model.rxns. The output from this
%   function does not include enzyme usage reactions, to observe these, on
%   could consider running flux variability directly on the ecModel.
%
% Input:
%   ecModel     an ecModel in GECKO 3 format (with ecModel.ec structure)
%   model       non-ecModel variant of the ecModel, to which the minFlux
%               and maxFlux will be mapped
% Output:
%   minFlux     vector of minimum flux rates, corresponding to model.rxns
%   maxFlux     vector of maximum flux rates, corresponding to model.rxns
%
% Usage: [minFlux, maxFlux] = ecFVA(ecModel, model)

rxnIDs = regexprep(ecModel.rxns,'(_REV)?(_EXP_\d+)?','');
[rxnIDmap, convRxnID] = findgroups(rxnIDs);

solMaxAll = nan(numel(ecModel.rxns),numel(convRxnID));
solMinAll = solMaxAll;

pool = gcp('nocreate');
if isempty(pool)
    parpool;
end

D = parallel.pool.DataQueue;
progressbar('Running ecFVA');
afterEach(D, @nUpdateProgressbar);

N = numel(convRxnID);
p = 1;

parfor i=1:N
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
        solMaxAll(:,i)=solMax.x;
    end
    tmpModel.c(rxnsToMax) = -1;
    tmpModel.c(rxnsToMin) = 1;
    solMin=solveLP(tmpModel);
    if ~isempty(solMin.x)
        solMinAll(:,i)=solMin.x;
    end    
    send(D, i);
end

minFlux=min(solMinAll,[],2,'omitnan');
maxFlux=max(solMaxAll,[],2,'omitnan');

mappedFlux = mapRxnsToConv(ecModel,model,[minFlux maxFlux]);

minFlux=mappedFlux(:,1);
maxFlux=mappedFlux(:,2);

function nUpdateProgressbar(~)
progressbar(p/N);
p = p + 1;
end
end
