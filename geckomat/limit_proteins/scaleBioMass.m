%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleBioMass(model,Ptot,GAM,scale_comp)
% 
% Benjamin Sanchez. Last update: 2018-10-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = scaleBioMass(model,Ptot,GAM,scale_comp)

if nargin < 3
    GAM = [];
end

%Option for changing composition & GAM (=true, default) or only GAM (=false):
if nargin < 4
    scale_comp = true;
end

%Compute carbohydrate and lipid new amounts, based on:
%1. Total mass remains constant, i.e. Pbase+Cbase+Lbase = Ptot+Ctot+Ltot
%2. Difference in mass is distributed proportionally, i.e. Ctot/Ltot = Cbase/Lbase
[~,Pbase,Cbase,~,~,Lbase] = sumBioMass(model);
Ctot = Cbase + (Pbase - Ptot)*Cbase/(Lbase+Cbase);
Ltot = Lbase + (Pbase - Ptot)*Lbase/(Lbase+Cbase);

%Compute rescaling fractions:
fP = Ptot/Pbase;
fC = Ctot/Cbase;
fL = Ltot/Lbase;

%Change compositions:
if scale_comp
    model = rescalePseudoReaction(model,'protein',fP);
    model = rescalePseudoReaction(model,'carbohydrate',fC);
    model = rescalePseudoReaction(model,'lipid backbone',fL);
    model = rescalePseudoReaction(model,'lipid chain',fL);
end

%Fit GAM if not available:
if isempty(GAM)
    GAM = fitGAM(model);
end

%Change GAM:
xr_pos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        %Polymerization costs from Forster et al 2003 - table S8:
        model.S(i,xr_pos) = sign(S_ix)*(GAM + 37.7*Ptot + 12.8*Ctot);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = rescalePseudoReaction(model,metName,f)

rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
for i = 1:length(model.mets)
    S_ir   = model.S(i,rxnPos);
    isProd = strcmp(model.metNames{i},metName);
    if S_ir ~= 0 && ~isProd
        model.S(i,rxnPos) = f*S_ir;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
