%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainPool(model,non_measured,UB)
% 
%
% Benjamín J. Sánchez. Last edited: 2016-04-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = constrainPool(model,non_measured,UB)

%Define new rxns: For each enzyme, add a new rxn that draws enzyme from the
%enzyme pool (a new metabolite), and remove previous exchange rxn. The new
%rxns have the following stoichiometry (T is the enzyme pool):
% MW[i]*P[T] -> P[i]
for i = 1:length(model.enzymes)
    if non_measured(i)
        rxnName = ['draw_prot_' model.enzymes{i}];
        metList = {'prot_pool' ['prot_' model.enzymes{i}]};
        model   = addReaction(model,rxnName,metList,[-model.MWs(i) 1],false,0,Inf);
        model   = removeRxns(model,['prot_' model.enzymes{i} '_exchange']);
    end
end

%Update metComps (last position is for the newly created protein pool):
model.metComps(end+1) = 2;      %For simplification the protein pool is in cytosol

%Finally, constraint enzyme pool by fixed value:
model = addReaction(model,'prot_pool_exchange',{'prot_pool'},1,false,0,UB);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%