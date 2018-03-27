%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeProtein(model,p,fs,options,GAM)
% 
%
% Benjamín J. Sánchez. Last edited: 2017-10-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeProtein(model,Ptot,fs,GAM,change_comp)

if nargin < 4
    GAM = 31;      %63.3% eff OXPHO - no H2O in prot/carb
end

%Option for changing composition & GAM (=true) or only GAM (=false):
if nargin < 5
    change_comp = true;
end

% Change Protein total amount:
Pbase = 0.4005;         %Value from biomass comp. (Förster data @ 0.1 1/h)
Cbase = 0.4067;         %Value from biomass comp. (Förster data @ 0.1 1/h)
P_pos = strcmp(model.rxns,'prot_pool_exchange');
if sum(P_pos) == 0
    disp('Metabolic model')
else
    model.ub(P_pos) = fs*Pbase;    %Fixed % protein
end

fP = Ptot/Pbase;
fC = (Cbase+Pbase-Ptot)/Cbase;

%Change biomass composition:
bio_pos = find(strcmp(model.rxns,'r_4041'));
for i = 1:length(model.mets)
    S_ix = model.S(i,bio_pos);
    if S_ix ~= 0        
        name  = model.metNames{i};
        isaa  = ~isempty(strfind(name,'tRNA'));
        isCH  = sum(strcmp({'(1->3)-beta-D-glucan','(1->6)-beta-D-glucan', ...
                            'chitin','glycogen','mannan','trehalose'},name)) == 1;
        isATP = strcmp(name,'ATP');
        isADP = strcmp(name,'ADP');
        isH2O = strcmp(name,'H2O');
        isH   = strcmp(name,'H+');
        isP   = strcmp(name,'phosphate');
        
        %Variable ATP growth related maintenance
        if isATP || isADP || isH2O || isH || isP
            S_ix = sign(S_ix)*(GAM + 16.965*fP + 5.210*fC);
        
        %Variable aa content in biomass eq:
        elseif isaa && change_comp
            S_ix = fP*S_ix;
            
        %Variable carb content in biomass eq:
        elseif isCH && change_comp
            S_ix = fC*S_ix;
        end
        
        model.S(i,bio_pos) = S_ix;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%