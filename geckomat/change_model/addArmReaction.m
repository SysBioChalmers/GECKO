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
model = addReaction(model,rxnID, ...
                    'reactionName', [model.rxnNames{rxnIndex} ' (arm)'], ...
                    'metaboliteList', [metS,['pmet_' rxn]], ...
                    'stoichCoeffList', [coeffsS,1], ...
                    'lowerBound', LB, ...
                    'upperBound', UB, ...
                    'objectiveCoef', obj, ...
                    'subSystem', subSystem);
model.grRules{strcmp(model.rxns,rxnID)} = grRule;

%Change old rxn:
model = addReaction(model,rxn, ...
                    'reactionName', model.rxnNames{rxnIndex}, ...
                    'metaboliteList', [['pmet_' rxn], metP], ...
                    'stoichCoeffList', [-1,coeffsP]);

%Update metComps:
pos = strcmp(model.mets,['pmet_' rxn]);
if sum(sub_pos) > 0
    model.metComps(pos) = model.metComps(sub_pos(1));
else
    model.metComps(pos) = model.metComps(pro_pos(1));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%