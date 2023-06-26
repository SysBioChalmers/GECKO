function [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel,proteinId)
% getReactionsFromEnzyme
%   Get all reactions that are annotated to a particular enzyme.
%
% Input:
%   ecModel     an ecModel in GECKO 3 format (with ecModel.ec structure)
%   proteinId   protein identifier, matching ecModel.ec.enzymes.
%
% Output:
%   rxns        reactions that are associated with this enzyme
%   kcat        kcat values of the corresponding reactions
%   idx         index of the reactions in ecModel.ec.rxns
%   rxnNames    names of the reactions
%   grRules     grRules of the reactions
%
% Usage: [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel,proteinId)

protIdx     = find(strcmpi(ecModel.ec.enzymes,proteinId));
ecRxnIdx    = find(ecModel.ec.rxnEnzMat(:,protIdx));
rxns        = ecModel.ec.rxns(ecRxnIdx);
kcat        = ecModel.ec.kcat(ecRxnIdx);

idx         = ecRxnIdx;

[~,rxnIdx]  = ismember(rxns,ecModel.rxns);
rxnNames    = ecModel.rxnNames(rxnIdx);
grRules     =  ecModel.grRules(rxnIdx);
end
