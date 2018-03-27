%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = modelCorrections(model)
% Corrects various potential issues in an EC-GEM
%
% INPUT:    A GEM model as a .mat structure.
% OUTPUT:   The corrected model.
%
% Ivan Domenzain. Last edited: 2018-02-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modelCorrectionsYli(model)
    %Enclose every isoenzyme into brackets
    for i=1:length(model.grRules)
        STR   = model.grRules{i};
        STR   = strrep(STR,' or ',' OR ');
        STR   = strrep(STR,' and ',' AND ');
        ORpos = strfind(STR,' OR ');
        if ~isempty(ORpos)
            for j=1:length(ORpos)
                ORpos = strfind(STR,' OR ');
                tempSTR = STR(ORpos(j):end);
                if ~strcmpi(STR(ORpos(j)-1),')')
                    STR = horzcat(STR(1:ORpos(j)-1),')',tempSTR);
                end
                ORpos   = strfind(STR,' OR ');
                tempSTR = STR(1:ORpos(j)+3);
                if ~strcmpi(STR(ORpos(j)+4),'(')
                    STR = horzcat(tempSTR,'(',STR(ORpos(j)+4:end));
                end
            end
        end
        model.grRules{i}  = STR;        
    end
    %Remove compartments from mets and create metComps field with this
    %information
    if isfield(model,'mets')
        if ~isfield(model,'metComps')
            model.metComps = ones(length(model.mets),1);
        end
        
        for i=1:length(model.mets)
            STR = model.mets{i};
            if strcmp(STR(end-2),'[') && strcmp(STR(end),']')
                comp  = STR(end-1);
                index = find(strcmpi(model.comps,comp),1);
                if ~isempty(index)
                    model.metComps(i) = index;
                end
                model.mets{i} = STR(1:end-3);
            end
            
        end
    end
                
            
end