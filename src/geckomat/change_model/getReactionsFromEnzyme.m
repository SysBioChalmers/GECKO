function [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel,proteinId)
% getReactionsFromEnzyme  Get all reactions annotated to a particular enzyme.
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% proteinId : char
%     protein identifier, matching ecModel.ec.enzymes case-insensitively.
%
% Returns
% -------
% rxns : cell
%     reactions that are associated with this enzyme.
% kcat : double
%     kcat values of the corresponding reactions.
% idx : double
%     index of the reactions in ecModel.ec.rxns.
% rxnNames : cell
%     names of the reactions.
% grRules : cell
%     grRules of the reactions.
%
% Notes
% -----
% An unmatched proteinId is not an error: all five outputs are returned
% empty, with no warning. Matching against ecModel.ec.enzymes is
% case-insensitive.
%
% Examples
% --------
%     [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel, proteinId);

protIdx     = find(strcmpi(ecModel.ec.enzymes,proteinId));
ecRxnIdx    = find(ecModel.ec.rxnEnzMat(:,protIdx));
rxns        = ecModel.ec.rxns(ecRxnIdx);
kcat        = ecModel.ec.kcat(ecRxnIdx);

idx         = ecRxnIdx;

[~,rxnIdx]  = ismember(rxns,ecModel.rxns);
rxnNames    = ecModel.rxnNames(rxnIdx);
grRules     =  ecModel.grRules(rxnIdx);
end
