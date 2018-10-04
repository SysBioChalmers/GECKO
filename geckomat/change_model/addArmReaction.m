%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addArmReaction(model,reaction)
% Adds an arm reaction for selected reaction in the model
%
% INPUT:
% model             The GEM structure (1x1 struct)
% reaction          The reaction ID   (string)
%
% OUTPUTS:
% model             Modified GEM structure (1x1 struct)
% 
% Cheng Zhang & Ivan Domenzain. Last edited: 2018-05-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addArmReaction(model,rxn)

%Define all neccesary parts for new rxn:
rxnIndex = strcmp(model.rxns,rxn);
sub_pos  = find(model.S(:,rxnIndex) < 0);
pro_pos  = find(model.S(:,rxnIndex) > 0);
metS     = model.mets(sub_pos)';
metP     = model.mets(pro_pos)';
LB       = model.lb(rxnIndex);
UB       = model.ub(rxnIndex);
obj      = model.c(rxnIndex);
coeffsS  = model.S(sub_pos,rxnIndex)';
coeffsP  = model.S(pro_pos,rxnIndex)';
grRule   = model.grRules{rxnIndex};

subSystem = '';
if isfield(model,'subSystems')
    if ~isempty(model.subSystems{rxnIndex}{1})
        subSystem = model.subSystems{rxnIndex};
    end
end

%Create new rxn:
rxnID = ['arm_' rxn];
rxnsToAdd.rxns = {rxnID};
rxnsToAdd.rxnNames = {[model.rxnNames{rxnIndex} ' (arm)']};
rxnsToAdd.mets = [metS,['pmet_' rxn]];
rxnsToAdd.stoichCoeffs = [coeffsS,1];
rxnsToAdd.lb = LB;
rxnsToAdd.ub = UB;
rxnsToAdd.obj = obj;
if isfield(model,'subSystems')
    rxnsToAdd.subSystems = subSystem;
end
model = addRxns(model,rxnsToAdd,1,'c',true); % All 'arm' metabolites are initially located to the cytosol
model.grRules{strcmp(model.rxns,rxnID)} = grRule;

%Change old rxn:
equations.mets = [['pmet_' rxn], metP];
equations.stoichCoeffs = [-1,coeffsP];
model = changeRxns(model,rxn,equations,1,'c');

%Update metComps:
pos = strcmp(model.mets,['pmet_' rxn]);
if sum(sub_pos) > 0
    model.metComps(pos) = model.metComps(sub_pos(1));
else
    model.metComps(pos) = model.metComps(pro_pos(1));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%