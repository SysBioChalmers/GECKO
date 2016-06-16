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
% Cheng Zhang. Last edited: 2016-03-18
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

%Create new rxn:
mets   = [metS,['pmet_' rxn]];
coeffs = [coeffsS,1];
model  = addReaction(model,['arm_' rxn],mets,coeffs,true,LB,UB,obj);

%Change old rxn:
model  = addReaction(model,rxn,[['pmet_' rxn], metP],[-1,coeffsP]);

%Update metComps (last position is for the newly created pmet):
model.metComps(end+1) = model.metComps(sub_pos(1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%