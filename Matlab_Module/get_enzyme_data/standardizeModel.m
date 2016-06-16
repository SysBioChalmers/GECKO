%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = standardizeModel(model,toolbox)
% Standardize a yeast GEM model to include some RAVEN fields.
%
% INPUT:    A yeast model as a .mat structure and the toolbox name 
%           used to create said model ('RAVEN' or 'COBRA').
% OUTPUT:   The standardized model. Main corrections applied are:
%           *Removal of any compartment reference in MetNames
%           *Addition of the comps vector (if not present)
%           *Addition of the compNames vector (if not present)
%           *Addition of the metComps vector (if not present)
%
% Benjamín J. Sánchez. Last edited: 2015-09-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = standardizeModel(model,toolbox)

mets = model.metNames;
M    = length(mets);

if strcmp(toolbox,'COBRA')
    comps = cell(M,1);
    for i = 1:M
        name     = mets{i};
        m        = length(name);
        %Split the met name in 2: The compartment and the metabolite name.
        %Afterwards, assign each to the corresponding list.
        k        = strfind(name,'[');
        n        = length(k);
        mets{i}  = name(1:k(n)-2);
        comps{i} = lower(name(k(n)+1:m-1));
    end
    model.metNames = mets;
    
    %Create compNames (list with the copartment names), comps
    %(abbreviations) and metComps (compartment in which each metabolite
    %is, using indexing of compNames):
    [comps_u,~] = deleteRepeated(comps);
    comps_c     = comps_u;
    for i = 1:length(comps_u)
        comp_name  = comps_u{i};
        comps_c{i} = comp_name(1);
        spaces = strfind(comp_name,' ');
        for j = spaces
            comps_c{i} = [comps_c{i} comp_name(j+1)];
        end
    end
        
    met_comps   = zeros(M,1);
    for i = 1:M
        for j = 1:length(comps_u)
            if strcmp(comps{i},comps_u{j})
                met_comps(i) = j;
            end
        end
    end
    model.comps     = comps_c;
    model.compNames = comps_u;
    model.metComps  = met_comps;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%