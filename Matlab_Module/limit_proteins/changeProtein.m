%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeProtein(model,p,fs,options,GAM)
% 
%
% Benjamín J. Sánchez. Last edited: 2016-03-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeProtein(model,P,fs,options,GAM)

%Current values (aerobic Yeast 7.6)
if nargin < 5
    GAM = 41.93;
    if nargin < 4
        options = [2;3;0.7;1;2;2;0.4266];
    end
end

% Change Protein total amount:
Pbase = options(7);
Cbase = 0.4067;
P_pos = strcmp(model.rxns,'prot_pool_exchange');
if sum(P_pos) == 0
    disp('Metabolic model')
else
    if options(4) == 1
        model.ub(P_pos) = fs*Pbase;    %Fixed % protein
    else
        model.ub(P_pos) = fs*P;        %Variable % protein
    end
end

Pfactor = P/Pbase;
Cfactor = (Cbase+Pbase-P)/Cbase;

%Change biomass composition:
bio_pos = find(strcmp(model.rxns,'r_4041'));
for i = 1:length(model.mets)
    S_ix = model.S(i,bio_pos);
    if S_ix ~= 0        
        name  = model.metNames{i};
        isaa  = ~isempty(strfind(name,'tRNA'));
        isATP = strcmp(name,'ATP');
        isADP = strcmp(name,'ADP');
        isH2O = strcmp(name,'H2O');
        isH   = strcmp(name,'H');
        isP   = strcmp(name,'phosphate');
        isCH  = ~isempty(strfind(name,'D-glucan')) || ...
                strcmp(name,'chitin') || strcmp(name,'glycogen') || ...
                strcmp(name,'mannan') || strcmp(name,'trehalose');
        
        %Variable ATP growth related maintenance
        if isATP || isADP || isH2O || isH || isP
            correction = 0.029*isH - 0.576*isP;
            if options(2) == 1                  %All constant
                S_ix = sign(S_ix)*(GAM + correction + 16.965 + 5.210);
            elseif options(2) == 2              %All is variable
                S_ix = sign(S_ix)*(GAM + correction + 16.965 + 5.210)*Pfactor;
                
            else                                %Only polim. is variable
                S_ix = sign(S_ix)*(GAM + correction + 16.965*Pfactor + 5.210*Cfactor);
            end
        
        %Variable aa content in biomass eq:
        elseif isaa
            if options(5) == 2
                S_ix = S_ix*Pfactor;
            end
            
        %Variable carb content in biomass eq:
        elseif isCH
            if options(6) == 2
                S_ix = S_ix*Cfactor;
            end
        end
        
        model.S(i,bio_pos) = S_ix;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%