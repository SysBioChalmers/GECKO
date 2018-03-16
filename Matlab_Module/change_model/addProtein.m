%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addProtein(model,P,kegg,swissprot)
% Adds an exchange reaction for protein P and updates model.enzymes,
% model.MWs and model.pathways to account for P.
%
% INPUT:
% model             Model with enzymes
% p                 Uniprot code of the protein
% kegg              KEGG database
% swissprot         Swissprot database
%
% OUTPUTS:
% model             Model with the added protein
% 
% Cheng Zhang & Benjamín Sánchez. Last edited: 2018-03-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addProtein(model,P,kegg,swissprot)

%Update model.enzyme vector:
prot_name     = ['prot_' P];
model.enzymes = [model.enzymes;P];
pos_e         = strcmp(model.enzymes,P);        %position in model.enzymes

%Update model.MWs & model.sequences vectors:
match_geneName = false;
match_MW       = false;
match_seq      = false;
for i = 1:length(swissprot)
    %Gene name:
    if strcmp(P,swissprot{i,1}) && ~isempty(swissprot{i,3}) && ~match_geneName
        match_geneName  = true;
        geneName        = swissprot{i,3};
        pos_space       = strfind(geneName,' ');
        if isempty(pos_space)
            pos_space = length(geneName)+1;
        end
        model.geneNames{pos_e,1} = geneName(1:pos_space(1)-1);
    end
    %Molecular Weight:
    if strcmp(P,swissprot{i,1}) && swissprot{i,5} ~= 0 && ~match_MW
        match_MW           = true;
        model.MWs(pos_e,1) = swissprot{i,5}/1000;	%g/mmol
    end
    %Sequence:
    if strcmp(P,swissprot{i,1}) && ~isempty(swissprot{i,6}) && ~match_seq
        match_seq                = true;
        model.sequences{pos_e,1} = swissprot{i,6};
    end
end
if ~match_geneName
    model.geneNames{pos_e,1} = '-';
end
if ~match_MW
    model.MWs(pos_e,1) = mean(cell2mat(swissprot(:,5)))/1000;	%average g/mmol
end
if ~match_seq
    model.sequences{pos_e,1} = '-';
end

%Update model.genes & model.pathways vectors:
match_gen  = false;
match_path = false;
for i = 1:length(kegg)
    if strcmp(P,kegg{i,1})
        %Gene:
        if ~isempty(kegg{i,3}) && ~match_gen
            match_gen            = true;
            model.genes{pos_e,1} = kegg{i,3};
        end
        %Pathway:
        if ~isempty(kegg{i,6}) && ~match_path
            match_path              = true;
            model.pathways{pos_e,1} = kegg{i,6};
        end
        %Molecular Weight (if nothing found in uniprot):
        if kegg{i,5} > 0 && ~match_MW
            match_MW           = true;
            model.MWs(pos_e,1) = kegg{i,5}/1000;
        end
        %Sequence (if nothing found in uniprot):
        if ~isempty(kegg{i,7}) && ~match_seq
            match_seq                = true;
            model.sequences(pos_e,1) = kegg{i,7};
        end
    end
end
if ~match_gen
    unknowns = ~cellfun(@isempty,strfind(model.genes,'unknown_'));
    if sum(unknowns) == 0
        idx = 0;
    else
        unknowns  = model.genes(unknowns);
        pos_final = strfind(unknowns{end},'_')+1;
        idx       = str2double(unknowns{end}(pos_final:end));
    end
    model.genes{pos_e,1} = ['unknown_' num2str(idx+1)];
end
if ~match_path
    model.pathways{pos_e,1} = '-';
end

%Add exchange reaction of protein: -> P
model = addReaction(model, ...                      %model
                    ['prot_' P '_exchange'], ...    %rxn name
                    {prot_name}, ...                %metabolites
                    1, ...                          %stoichiometry
                    true, ...                       %reversibility
                    0, ...                          %LB
                    Inf, ...                        %UB
                    0, ...                          %c
                    {''}, ...                       %subsystem
                    model.genes{pos_e,1});          %gene rule

%Update metComps:
pos_m = strcmp(model.mets,prot_name);   %position in model.mets
if isfield(model,'compNames')
    cytIndex = find(strcmpi(model.compNames,'cytoplasm'),1);
    if ~isempty(cytIndex)
        model.metComps(pos_m) = 2;              %For simplification all proteins are in cytosol
    else
        model.metComps(pos_m) = 1; 
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%