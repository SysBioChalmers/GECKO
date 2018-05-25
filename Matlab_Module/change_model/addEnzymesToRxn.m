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
%
% OUTPUTS:
% model             Modified GEM structure (1x1 struct)
% 
% Cheng Zhang & Ivan Domenzain. Last edited: 2018-05-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addEnzymesToRxn(model,kvalues,rxn,newMets,newRxnName,kegg,swissprot)

%Define all neccesary parts for new (or changed) rxn:
rxnIndex  = find(ismember(model.rxns,rxn)); 
metS      = model.mets(model.S(:,rxnIndex) < 0)';
metP      = model.mets(model.S(:,rxnIndex) > 0)';
LB        = model.lb(rxnIndex);
UB        = model.ub(rxnIndex);    
obj       = model.c(rxnIndex);
coeffsS   = model.S(model.S(:,rxnIndex)<0,rxnIndex)';
coeffsP   = model.S(model.S(:,rxnIndex)>0,rxnIndex)';
genes     = cell(size(newMets));

if isfield(model,'subSystems')
    subSystem = model.subSystems(rxnIndex);
else
    subSystem = '';
end

%Find genes either in swissprot or in kegg and with them construct the gene rule:
for i = 1:length(newMets)
    protein = strrep(newMets{i},'prot_','');
    try
        geneIDs      = swissprot{strcmp(swissprot(:,1),protein),3};
        geneIDs      = strsplit(geneIDs,' ');
        [genes(i),~] = intersect(geneIDs,model.genes);
    catch
        genes{i} = kegg{strcmp(kegg(:,1),protein),3};
    end
end
grRule = strjoin(genes,' and ');

%Include enzyme in reaction:
mets   = [metS,newMets,metP];
coeffs = [coeffsS,-kvalues.^-1,coeffsP];
model  = addReaction(model,newRxnName,mets,coeffs,true,LB,UB,obj,subSystem,grRule,'','',false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%