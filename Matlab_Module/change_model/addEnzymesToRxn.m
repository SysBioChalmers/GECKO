%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName)
% Adds new metabolite to the left side of a selected reaction in the model.
% If the reaction does not exist it will create a new one.
%
% INPUT:
% model        the GEM structure (1x1 struct)
% kvalues      kcat values of the enzyme/complex
% rxn          the reaction original ID   (string)
% newMets      name of the new pseudo-metabolite (enzyme)
% newRxnName   {ID,name} of the new reaction
%
% OUTPUTS:
% model             Modified GEM structure (1x1 struct)
% 
% Cheng Zhang. Last edited: 2016-03-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName,n)

%Define all neccesary parts for new (or changed) rxn:
rxnIndex  = find(ismember(model.rxns,rxn)); 
metS      = model.mets(model.S(:,rxnIndex) < 0)';
metP      = model.mets(model.S(:,rxnIndex) > 0)';
LB        = model.lb(rxnIndex);
UB        = model.ub(rxnIndex);    
obj       = model.c(rxnIndex);
coeffsS   = model.S(model.S(:,rxnIndex)<0,rxnIndex)';
coeffsP   = model.S(model.S(:,rxnIndex)>0,rxnIndex)';
grRule    = model.grRules(rxnIndex);
if isfield(model,'subSystems')
    subSystem = model.subSystems(rxnIndex);
else
    subSystem = '';
end

%Include enzyme in reaction:
for i = 1:length(newMets)
    metS    = [metS,newMets{i}];
    coeffsS = [coeffsS,-1/kvalues(i)];
end
mets   = [metS,metP];
coeffs = [coeffsS,coeffsP];
model  = addReaction(model,newRxnName,mets,coeffs,true,LB,UB,obj,subSystem,'','','',false);
if n==1
    %For reactions with no isoenzymes
    model.grRules(end) = grRule;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%