function [minFlux, maxFlux] = ecFVA(ecModel, model)
% ecFVA
%   Flux variability analysis is performed on the ecModel, and isoenzymic
%   reactions are combined to construct ouput minFlux and maxFlux vectors,
%   which follow the same order of model.rxns. The output from this
%   function does not include enzyme usage reactions, to observe these, on
%   could consider running flux variability directly on the ecModel.
%
%    ecModel    will be used in the flux variability analysis
%    model      non-ecModel variant of the ecModel, to which the minFlux
%               and maxFlux will be mapped
%
%    minFlux    vector of minimum flux rates, corresponding to model.rxns
%    maxFlux    vector of maximum flux rates, corresponding to model.rxns
%
% Usage: [minFlux, maxFlux] = ecFVA(ecModel, model)

rxnIDs = regexprep(ecModel.rxns,'(_REV)?(_EXP_\d+)?','');
[rxnIDmap, convRxnID] = findgroups(rxnIDs);

solMaxAll = zeros(numel(ecModel.rxns),numel(convRxnID));
solMinAll = solMaxAll;

D = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...');
afterEach(D, @nUpdateWaitbar);

N = numel(rxnIDmap);
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
    solMaxAll(:,i)=solMax.x;
    tmpModel.c(rxnsToMax) = -1;
    tmpModel.c(rxnsToMin) = 1;
    solMin=solveLP(tmpModel);
    solMinAll(:,i)=solMin.x;
    send(D, i);
end

minFlux=min(solMinAll,[],2);
maxFlux=max(solMaxAll,[],2);

mappedFlux = mapRxnsToConv(ecModel,model,[minFlux maxFlux]);

minFlux=mappedFlux(:,1);
maxFlux=mappedFlux(:,2);
function nUpdateWaitbar(~)
waitbar(p/N, h);
p = p + 1;
end
end


