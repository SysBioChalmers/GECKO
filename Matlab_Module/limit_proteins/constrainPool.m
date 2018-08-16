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
        rxnID = ['draw_prot_' model.enzymes{i}];
        model = addReaction(model,rxnID, ...
                            'metaboliteList', {'prot_pool' ['prot_' model.enzymes{i}]}, ...
                            'stoichCoeffList', [-model.MWs(i) 1], ...
                            'reversible', false, ...
                            'lowerBound', 0, ...
                            'upperBound', Inf);
        model.grRules{strcmp(model.rxns,rxnID)} = model.enzGenes{i};
        model = removeReactions(model,{['prot_' model.enzymes{i} '_exchange']});
    end
end

%Finally, constraint enzyme pool by fixed value:
model = addReaction(model,'prot_pool_exchange', ...
                    'metaboliteList', {'prot_pool'}, ...
                    'stoichCoeffList', 1, ...
                    'reversible', false, ...
                    'lowerBound', 0, ...
                    'upperBound', UB);

%Update metComps (last position is for the newly created protein pool):
cytIndex = find(strcmpi(model.compNames,'cytoplasm'),1);
if ~isempty(cytIndex)
    model.metComps(end) = cytIndex;	%For simplification all proteins are in cytosol
else
    model.metComps(end) = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%