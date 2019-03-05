function mappedIndxs = rxnMapping(rxnID,model,revFlag)
% rxnMapping
% 
% Function that maps a metabolic rxn from a GEM on its enzyme constrained 
% version (irreversible model).
%
%    rxnID         model.rxns(i) that is going to be mapped
%    model         ecModel in which the rxn is going to be searched
%    revFlag       True if the searched rxn is reversible
%
%    mappedIndxs   Cell array that contains the index of the corresponding
%                  metabolic rxn in the ecModel (the arm rxn in the case of
%                  isoenzymes). If the original rxn is reversible then it
%                  contains the indexes of the forward and backward rxns
%                  respectively.
%
% Usage: mappedIndxs = rxnMapping(rxnID,model,revFlag)
%
% Ivan Domenzain.      Last edited: 2019-03-04


indexes = find(contains(model.rxns,rxnID));
if revFlag
    backwardIndxs = indexes(find(contains(model.rxns(indexes),'_REV')));
    forwardIndxs  = setdiff(indexes,backwardIndxs);
    backwardIndxs = findArmRxns(backwardIndxs,model);
    forwardIndxs  = findArmRxns(forwardIndxs,model);
    mappedIndxs   = vertcat(forwardIndxs,backwardIndxs);
else
    mappedIndxs  = findArmRxns(indexes,model);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ArmIndex = findArmRxns(rxnIndexes,model)
if length(rxnIndexes)>1
    ArmIndex = rxnIndexes(contains(model.rxns(rxnIndexes),'arm_'));
else
    ArmIndex = rxnIndexes;
end
end