%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = deleteProtein(model,P)
%
% INPUT:
% model             Model with enzymes
% p                 Uniprot code of the protein
%
% OUTPUTS:
% model             Model with the deleted protein
% 
% Benjamín Sánchez. Last edited: 2018-05-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = deleteProtein(model,P)

%Find position:
enz_pos   = strcmp(model.enzymes,P);
met_pos   = strcmp(model.mets,['prot_' P]);
gene_name = model.enzGenes(enz_pos);

%Delete enzyme fields:
model.enzymes(enz_pos)   = [];
model.enzGenes(enz_pos)  = [];
model.enzNames(enz_pos)  = [];
model.MWs(enz_pos)       = [];
model.sequences(enz_pos) = [];
model.pathways(enz_pos)  = [];

%Delete metabolite also:
model = removeMets(model,model.mets(met_pos),true,true);

%Delete gene also:
model = removeGenes(model,gene_name);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%