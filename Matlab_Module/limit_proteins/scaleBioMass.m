%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleBioMass(model,Ptot,GAM,scale_comp)
% 
% Benjamín Sánchez. Last update: 2018-08-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = scaleBioMass(model,Ptot,GAM,scale_comp)

%Fix GAM if not provided:
if nargin < 3
    GAM = 31;      %63.3% eff OXPHO - no H2O in prot/carb
end

%Option for changing composition & GAM (=true, default) or only GAM (=false):
if nargin < 4
    scale_comp = true;
end

% Compute protein and carbohydrate fractions:
[~,Pbase,Cbase,~,~,~] = sumBioMass(model);
fP   = Ptot/Pbase;
Ctot = Cbase + Pbase - Ptot;
fC   = Ctot/Cbase;

if scale_comp
    %Change protein composition:
    pr_pos = strcmp(model.rxnNames,'protein pseudoreaction');
    for i = 1:length(model.mets)
        S_ip   = model.S(i,pr_pos);
        isProd = strcmp(model.metNames{i},'protein');
        if S_ip ~= 0 && ~isProd
            model.S(i,pr_pos) = fP*S_ip;
        end
    end
    
    %Change carbohydrate composition:
    cr_pos = strcmp(model.rxnNames,'carbohydrate pseudoreaction');
    for i = 1:length(model.mets)
        S_ic   = model.S(i,cr_pos);
        isProd = strcmp(model.metNames{i},'carbohydrate');
        if S_ic ~= 0 && ~isProd
            model.S(i,cr_pos) = fC*S_ic;
        end
    end
end

%Change GAM:
xr_pos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        %Polymerization costs from Förster et al 2003 - table S8:
        model.S(i,xr_pos) = sign(S_ix)*(GAM + 37.7*Ptot + 12.8*Ctot);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
