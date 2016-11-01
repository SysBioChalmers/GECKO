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
% Cheng Zhang. Last edited: 2016-10-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addProtein(model,P,kegg,swissprot)

%Add exchange reaction of protein: -> P
prot_name    = ['prot_' P];
exchange_rxn = ['prot_' P '_exchange'];
model        = addReaction(model,exchange_rxn,{prot_name},1,true,0,1000,0);

%Update model.enzyme vector:
model.enzymes = [model.enzymes;P];

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
        model.geneNames = [model.geneNames;geneName(1:pos_space(1)-1)];
    end
    %Molecular Weight:
    if strcmp(P,swissprot{i,1}) && swissprot{i,5} ~= 0 && ~match_MW
        match_MW  = true;
        model.MWs = [model.MWs;swissprot{i,5}/1000];	%g/mmol
    end
    %Sequence:
    if strcmp(P,swissprot{i,1}) && ~isempty(swissprot{i,6}) && ~match_seq
        match_seq      = true;
        model.sequence = [model.sequence;swissprot{i,6}];
    end
end
if ~match_geneName
    model.geneNames = [model.geneNames;'-'];
end
if ~match_MW
    model.MWs = [model.MWs;mean(cell2mat(swissprot(:,5)))/1000];	%average g/mmol
end
if ~match_seq
    model.sequence = [model.sequence;'-'];
end

%Update model.genes & model.pathways vectors:
match_gen  = false;
match_path = false;
for i = 1:length(kegg)
    if strcmp(P,kegg{i,1})
        %Gene:
        if ~isempty(kegg{i,3}) && ~match_gen
            match_gen    = true;
            model.genes2 = [model.genes2;kegg{i,3}];
        end
        %Pathway:
        if ~isempty(kegg{i,6}) && ~match_path
            match_path     = true;
            model.pathways = [model.pathways;kegg{i,6}];
        end
        %Molecular Weight (if nothing found in uniprot):
        if kegg{i,5} == 0 && ~match_MW
            match_MW       = true;
            model.MWs(end) = kegg{i,5}/1000;
        end
        %Sequence (if nothing found in uniprot):
        if ~isempty(kegg{i,7}) && ~match_seq
            match_seq           = true;
            model.sequence(end) = kegg{i,7};
        end
    end
end
if ~match_gen
    model.genes2 = [model.genes2;'-'];
end
if ~match_path
    model.pathways = [model.pathways;'-'];
end

%Update metComps (last position is for the newly created protein):
model.metComps(end+1) = 2;      %For simplification all proteins are in cytosol

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%