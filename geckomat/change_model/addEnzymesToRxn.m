%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName,protGenes)
% Adds new metabolite to the left side of a selected reaction in the model.
% If the reaction does not exist it will create a new one.
%
% INPUT:
% model        the GEM structure (1x1 struct)
% kvalues      kcat values of the enzyme/complex
% rxn          the reaction original ID   (string)
% newMets      name of the new pseudo-metabolite (enzyme)
% newRxnName   {ID,name} of the new reaction
% protGenes    Gene ID associated with the provided enzyme
%
% OUTPUTS:
% model             Modified GEM structure (1x1 struct)
% 
% Cheng Zhang & Ivan Domenzain. Last edited: 2018-09-07
% Eduard Kerkhoven & Benjamin Sanchez. Last edited: 2018-11-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName,protGenes)

if nargin < 5
    protGenes = '';
end

%Define all necessary parts for new (or changed) rxn:
rxnIndex = strcmp(model.rxns,rxn); 
metS     = model.mets(model.S(:,rxnIndex) < 0)';
metP     = model.mets(model.S(:,rxnIndex) > 0)';
coeffsS  = model.S(model.S(:,rxnIndex)<0,rxnIndex)';
coeffsP  = model.S(model.S(:,rxnIndex)>0,rxnIndex)';

%Find default compartment:
cytIndex = strcmpi(model.compNames,'cytoplasm');
if sum(cytIndex) == 1
    comp = model.comps{cytIndex};	%For simplification all proteins are in cytosol
else
    comp = model.comps{1};
end

%Include enzyme in reaction:
rxnToAdd.mets         = [metS,newMets,metP];
rxnToAdd.stoichCoeffs = [coeffsS,-kvalues.^-1,coeffsP];
if ismember(newRxnName{1},model.rxns)
    model = changeRxns(model,newRxnName(1),rxnToAdd,1,comp);
else    
    rxnToAdd.rxns     = newRxnName(1);
    rxnToAdd.rxnNames = newRxnName(2);
    rxnToAdd.lb       = model.lb(rxnIndex);
    rxnToAdd.ub       = model.ub(rxnIndex);
    rxnToAdd.obj      = model.c(rxnIndex);
    if ~isempty(protGenes)
        rxnToAdd.grRules = {protGenes};
    end
    if isfield(model,'subSystems')
        rxnToAdd.subSystems = model.subSystems(rxnIndex);
    end
    model = addRxns(model,rxnToAdd,1,comp,true);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
