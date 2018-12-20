%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainPool(model,non_measured,UB)
% 
% Benjamin J. Sanchez. Last edited: 2018-11-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = constrainPool(model,non_measured,UB)

%Find default compartment:
cytIndex = strcmpi(model.compNames,'cytoplasm');
if sum(cytIndex) == 1
    comp = model.comps{cytIndex};	%Protein pool in cytosol
else
    comp = model.comps{1};
end

%Define new rxns: For each enzyme, add a new rxn that draws enzyme from the
%enzyme pool (a new metabolite), and remove previous exchange rxn. The new
%rxns have the following stoichiometry (T is the enzyme pool):
% MW[i]*P[T] -> P[i]
for i = 1:length(model.enzymes)
    if non_measured(i)
        rxnToAdd.rxns         = {['draw_prot_' model.enzymes{i}]};
        rxnToAdd.rxnNames     = rxnToAdd.rxns;
        rxnToAdd.mets         = {'prot_pool' ['prot_' model.enzymes{i}]};
        rxnToAdd.stoichCoeffs = [-model.MWs(i) 1];
        rxnToAdd.lb           = 0; % ub is taken from model default, otherwise inf
        rxnToAdd.grRules      = model.enzGenes(i);
        model = addRxns(model,rxnToAdd,1,comp,true);
        model = removeReactions(model,{['prot_' model.enzymes{i} '_exchange']});
    end
end

%Finally, constraint enzyme pool by fixed value:
rxnToAdd.rxns         = {'prot_pool_exchange'};
rxnToAdd.rxnNames     = rxnToAdd.rxns;
rxnToAdd.mets         = {'prot_pool'};
rxnToAdd.stoichCoeffs = 1;
rxnToAdd.lb           = 0;
rxnToAdd.ub           = UB;
rxnToAdd.grRules      = {''};
model = addRxns(model,rxnToAdd);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
