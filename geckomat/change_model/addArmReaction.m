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
% Eduard Kerkhoven & Benjamin Sanchez. Last edited: 2018-11-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addArmReaction(model,rxn)

%Define all neccesary parts for new rxn:
rxnIndex = strcmp(model.rxns,rxn);
sub_pos  = find(model.S(:,rxnIndex) < 0);
pro_pos  = find(model.S(:,rxnIndex) > 0);
metS     = model.mets(sub_pos)';
metP     = model.mets(pro_pos)';
coeffsS  = model.S(sub_pos,rxnIndex)';
coeffsP  = model.S(pro_pos,rxnIndex)';

%Find default compartment:
if sum(sub_pos) > 0
    comp = model.comps{model.metComps(sub_pos(1))};
else
    comp = model.comps{model.metComps(pro_pos(1))};
end

%Create new rxn:
rxnToAdd.rxns         = {['arm_' rxn]};
rxnToAdd.rxnNames     = {[model.rxnNames{rxnIndex} ' (arm)']};
rxnToAdd.mets         = [metS,['pmet_' rxn]];
rxnToAdd.stoichCoeffs = [coeffsS,1];
rxnToAdd.lb           = model.lb(rxnIndex);
rxnToAdd.ub           = model.ub(rxnIndex);
rxnToAdd.obj          = model.c(rxnIndex);
rxnToAdd.grRules      = model.grRules(rxnIndex);
if isfield(model,'subSystems')
    rxnToAdd.subSystems = model.subSystems(rxnIndex);
end
model = addRxns(model,rxnToAdd,1,comp,true);

%Change old rxn:
equation.mets = [['pmet_' rxn], metP];
equation.stoichCoeffs = [-1,coeffsP];
model = changeRxns(model,rxn,equation,1,comp);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
