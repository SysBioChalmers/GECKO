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
% Cheng Zhang. Last edited: 2017-10-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = deleteProtein(model,P)

%Find position:
enz_pos = strcmp(model.enzymes,P);
met_pos = strcmp(model.mets,['prot_' P]);

%Delete enzyme fields:
model.enzymes(enz_pos)   = [];
model.genes(enz_pos)     = [];
model.geneNames(enz_pos) = [];
model.MWs(enz_pos)       = [];
model.sequences(enz_pos) = [];
model.pathways(enz_pos)  = [];

%Delete metabolite also:
model = removeMetabolites(model,model.mets(met_pos));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%