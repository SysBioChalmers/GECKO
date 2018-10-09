%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainPool(model,non_measured,UB)
% 
% Benjamín J. Sánchez. Last edited: 2018-08-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = constrainPool(model,non_measured,UB)

%Define new rxns: For each enzyme, add a new rxn that draws enzyme from the
%enzyme pool (a new metabolite), and remove previous exchange rxn. The new
%rxns have the following stoichiometry (T is the enzyme pool):
% MW[i]*P[T] -> P[i]
for i = 1:length(model.enzymes)
    if non_measured(i)
        rxnsToAdd.rxns          = {['draw_prot_' model.enzymes{i}]};
        rxnsToAdd.rxnNames      = rxnsToAdd.rxns;
        rxnsToAdd.mets          = {'prot_pool' ['prot_' model.enzymes{i}]};
        rxnsToAdd.stoichCoeffs  = [-model.MWs(i) 1];
        rxnsToAdd.lb            = 0; % ub is taken from model's default, otherwise inf
        model = addRxns(model,rxnsToAdd);
        model.grRules{strcmp(model.rxns,rxnID)} = model.enzGenes{i};
        model = removeReactions(model,{['prot_' model.enzymes{i} '_exchange']});
    end
end

%Finally, constraint enzyme pool by fixed value:
rxnsToAdd.rxns          = {'prot_pool_exchange'};
rxnsToAdd.rxnNames      = rxnsToAdd.rxns;
rxnsToAdd.mets          = {'prot_pool'};
rxnsToAdd.stoichCoeffs  = 1;
rxnsToAdd.lb            = 0;
rxnsToAdd.ub            = UB;
model = addRxns(model,rxnsToAdd);

%Update metComps (last position is for the newly created protein pool):
cytIndex = find(strcmpi(model.compNames,'cytoplasm'),1);
if ~isempty(cytIndex)
    model.metComps(end) = cytIndex;	%For simplification all proteins are in cytosol
else
    model.metComps(end) = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%