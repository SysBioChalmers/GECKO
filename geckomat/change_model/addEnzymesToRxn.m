%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName,kegg,swissprot)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName,protGenes)

%Define all neccesary parts for new (or changed) rxn:
rxnIndex = strcmp(model.rxns,rxn); 
metS     = model.mets(model.S(:,rxnIndex) < 0)';
metP     = model.mets(model.S(:,rxnIndex) > 0)';
LB       = model.lb(rxnIndex);
UB       = model.ub(rxnIndex);    
obj      = model.c(rxnIndex);
coeffsS  = model.S(model.S(:,rxnIndex)<0,rxnIndex)';
coeffsP  = model.S(model.S(:,rxnIndex)>0,rxnIndex)';

subSystem = {''};
if isfield(model,'subSystems')
    if ~isempty(model.subSystems{rxnIndex}{1})
        subSystem = model.subSystems(rxnIndex);
    end
end

%Include enzyme in reaction:
rxnsToAdd.rxns = newRxnName(1);
rxnsToAdd.rxnNames = newRxnName(2);
rxnsToAdd.mets = [metS,newMets,metP];
rxnsToAdd.stoichCoeffs = [coeffsS,-kvalues.^-1,coeffsP];
rxnsToAdd.lb = LB;
rxnsToAdd.ub = UB;
rxnsToAdd.obj = obj;
if isfield(model,'subSystems')
    rxnsToAdd.subSystems = subSystem;
end
model = addRxns(model,rxnsToAdd,1,'c',true);

%Add/modify gene(s) association information if available
if nargin >5 && ~isempty(protGenes)
    model.grRules{strcmp(model.rxns,newRxnName{1})} = protGenes;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%