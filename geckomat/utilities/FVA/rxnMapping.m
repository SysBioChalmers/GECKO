%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function mappedIndxs = rxnMapping(rxnID,model,revFlag)
%  
% Function that maps a metabolic rxn from a GEM on its enzyme constrained 
% version (irreversible model).
%
% INPUTS:
%   - rxnID         model.rxns(i) that is going to be mapped
%   - model         ecModel in which the rxn is going to be searched
%   - revFlag       True if the searched rxn is reversible
% OUTPUTS:
%   - mappedIndxs   Cell array that contains the index of the corresponding
%                   metabolic rxn in the ecModel (the arm rxn in the case of
%                   isoenzymes). If the original rxn is reversible then it
%                   contains the indexes of the forward and backward rxns
%                   respectively.
%
% Ivan Domenzain.      Last edited: 2018-03-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mappedIndxs = rxnMapping(rxnID,model,revFlag)
    indexes = find(contains(model.rxns,rxnID));
    if revFlag
        backwardIndxs = find(contains(model.rxns(indexes),'_REV'));
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