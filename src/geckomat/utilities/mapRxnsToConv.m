function [mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, model, fluxVect)
% mapRxnsToConv
%   A vector (or matrix) of fluxes is mapped to the reactions in the
%   conventional starting model that was used to construct the ecModel. It
%   is essential that the provided conventional model is indeed the model
%   that was used to initate ecModel reconstruction.
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure),
%                   that was used to obtain fluxVect
%   model           the starting model for ecModel, to which the reactions
%                   should be mapped
%   fluxVect        vector or matrix of flux values, matching ecModel.rxns
%
% Output:
%   mappedFlux      vector or matrix of flux values, matching model.rxns
%   enzUsageFlux    vector or matrix of flux values from enzyme usage
%                   reactions, as these are absent from mappedFlux
%   usageEnz        cell array with protein IDs, matching enzUsageFlux
%
% Usage:
%   [mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, model, fluxVect)

fluxes = fluxVect;
rxnIDs = ecModel.rxns;

% Invert flux of _REV reactions
revRxns = endsWith(rxnIDs,'_REV') | contains(rxnIDs,'_REV_EXP_');
fluxes(revRxns,:) = -fluxes(revRxns,:);
rxnIDs(revRxns) = replace(rxnIDs(revRxns),'_REV','');
% Remove _EXP_. suffixes
rxnIDs = regexprep(rxnIDs,'_EXP_\d+','');

% Map and sum fluxes to converted reaction IDs
[rxnIDmap, convRxnID] = findgroups(rxnIDs);
newVect = splitapply(@(x){sum(x,1)}, fluxes, rxnIDmap);
newVect = cell2mat(newVect);

% Place in same order as in original model
[mapCheck,origIdx] = ismember(model.rxns,convRxnID);
if ~all(mapCheck)
    error('Not all reactions from model.rxns can be found in the ecModel. Are you sure that ecModel is derived from model?')
end
mappedFlux=newVect(origIdx,:);

% Separately report enzyme usages
usageEnz = startsWith(ecModel.rxns,{'usage_prot_','prot_pool_exchange'});
enzUsageFlux = fluxVect(usageEnz,:);
usageEnz = regexprep(ecModel.rxns(usageEnz),'usage_prot_','');
usageEnz = regexprep(usageEnz,'prot_pool_exchange','pool');
end
