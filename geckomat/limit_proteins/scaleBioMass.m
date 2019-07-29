%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleBioMass(model,Ptot,GAM,scale_comp)
% 
% Benjamin Sanchez. Last update: 2018-10-23
% Ivan Domenzain.   Last update: 2019-07-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = scaleBioMass(model,Ptot,GAM,scale_comp)
if nargin < 3
    GAM = [];
end
%Option for changing composition & GAM (=true, default) or only GAM (=false):
if nargin < 4
    scale_comp = true;
end
%Get biomass pseudoreaction ID and biomass components pseudoreactions names
cd ..
parameters = getModelParameters;
cd limit_proteins
%Compute carbohydrate and lipid new amounts, based on:
%1. Total mass remains constant, i.e. Pbase+Cbase+Lbase = Ptot+Ctot+Ltot
%2. Difference in mass is distributed proportionally, i.e. Ctot/Ltot = Cbase/Lbase
[~,Pbase,Cbase,R,D,Lbase] = sumBioMass(model);
Ctot = Cbase + (Pbase - Ptot)*Cbase/(Lbase+Cbase);
Ltot = Lbase + (Pbase - Ptot)*Lbase/(Lbase+Cbase);
%Compute rescaling fractions:
fP = Ptot/Pbase;
fC = Ctot/Cbase;
fL = Ltot/Lbase;
%Change compositions:
if scale_comp
    if isfield(parameters,'bio_comp')
        model = rescalePseudoReaction(model,parameters.bio_comp{1},fP);
        model = rescalePseudoReaction(model,parameters.bio_comp{2},fC);
        model = rescalePseudoReaction(model,parameters.bio_comp{3},fL);
        model = rescalePseudoReaction(model,parameters.bio_comp{4},fL);
    end
end
%Fit GAM if not available:
if isempty(GAM)
    GAM = fitGAM(model);
end
%Change GAM:
xr_pos = strcmp(model.rxns,parameters.bioRxn);
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        if strcmpi(parameters.keggID,'sce')
            %Polymerization costs from Forster et al 2003 - table S8:
            GAMpol = Ptot*37.7 + Ctot*12.8 + R*26.0 + D*26.0;
        end
        model.S(i,xr_pos) = sign(S_ix)*(GAM + GAMpol);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = rescalePseudoReaction(model,metName,f)
rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
if sum(rxnPos) == 1
    for i = 1:length(model.mets)
        S_ir   = model.S(i,rxnPos);
        isProd = strcmp(model.metNames{i},metName);
        if S_ir ~= 0 && ~isProd
            model.S(i,rxnPos) = f*S_ir;
        end
    end
end
end