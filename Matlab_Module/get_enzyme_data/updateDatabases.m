%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateDatabases
% Updates all databases for protein matching (KEGG and Swiss-Prot).
%
% Note: Before using this script, one should manually download:
%       *KEGG:      
%       *Swissprot: Download a tab delimited file with the following format:
%                   Entry - Protein names - Gene names - EC number - Sequence
%                   http://www.uniprot.org/uniprot/?query=organism:%22yeast%22
%                   OBS: filter with the Swiss-Prot option
% 
% Benjamín Sánchez. Last edited: 2016-01-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateDatabases

%Retrieve Swissprot Data (Entry - Protein names - Gene names - EC number - Sequence):
cd ../../Databases
fileID_uni        = fopen('uniprot-organism%3Ayeast.tab');
swissprot         = textscan(fileID_uni,'%s %s %s %s %s %s','delimiter','\t');
swissprot         = [swissprot{1} swissprot{2} swissprot{3} swissprot{4} swissprot{5}];
swissprot(1,:)    = [];
swissprot_complex = cell(1000,3);
fclose(fileID_uni);
cd ../Matlab_Module/get_enzyme_data
m = 0;
for i = 1:length(swissprot)
    %Leave protein name as lower case, remove ';' from ECs & calculate MW:
    prot_name      = lower(swissprot{i,2});
    uni            = swissprot{i,1};
    sequence       = swissprot{i,5};
    MW             = calculateMW(sequence);
    swissprot{i,2} = prot_name;
    swissprot{i,4} = strrep(swissprot{i,4},';','');
    swissprot{i,5} = MW;
    swissprot{i,6} = sequence;
    %Update complex database:
    [swissprot_complex,m] = updateComplexDB(swissprot_complex,m,prot_name,uni,MW);
    disp(['Updating Swiss-Prot database: Ready with protein ' uni])
end
swissprot_complex(m+1:end,:) = [];

%Retrieve KEGG info (gen - uniprot - EC code - name - MW - Pathway):
cd ../../Databases/KEGG
file_names      = dir();
file_names(1:2) = [];
kegg            = cell(100000,7);
kegg_complex    = cell(1000,3);
n               = 0;
m               = 0;
for i = 1:length(file_names)
    file_name = file_names(i).name;
    %1st column: Gene name
    gene_name = file_name(1:end-4);
    %Retrieve all data as a cell with all rows:
    fID  = fopen(file_name);
    text = textscan(fID,'%s','delimiter','\t');
    fclose(fID);
    text = text{1};
    cd ../../Matlab_Module/get_enzyme_data
    
    uni      = '';
    sequence = '';
    MW       = 0;
    pathway  = '';
    for j = 1:length(text)
        line = text{j};
        %2nd column: uniprot number
        if ~isempty(strfind(line,'UniProt:'))
            uni = line(10:end);
            
        %3rd & 4th column: protein name and EC number
        elseif ~isempty(strfind(line,'DEFINITION'))
            pos_EC    = strfind(line,'EC:');
            if isempty(pos_EC)
                prot_name = lower(line(13:end));
                EC_names  = '';
            else
                prot_name = lower(line(13:pos_EC-3));
                EC_names  = line(pos_EC+3:end-1);
            end
            
        %5th column and 7th column: MW & sequence
        elseif ~isempty(strfind(line,'AASEQ'))
            end_seq  = false;
            for k = j+1:length(text)
                if ~isempty(strfind(text{k},'NTSEQ'))
                    end_seq = true;
                elseif ~end_seq
                    sequence = [sequence text{k}];
                end
            end
            MW = calculateMW(sequence);
        
        %6th column: Pathway
        elseif ~isempty(strfind(line,'PATHWAY'))
            start    = strfind(line,'sce');
            pathway  = line(start(1):end);
            end_path = false;
            for k = j+1:length(text)
                nospace = strrep(text{k},'sce01100  Metabolic pathways','');
                nospace = strrep(nospace,' ','');
                if length(nospace) > 10
                    if strcmp(nospace(1:3),'sce') && ~end_path
                        start    = strfind(text{k},'sce');
                        pathway  = [pathway ' ' text{k}(start(1):end)];
                    else
                        end_path = true;
                    end
                end
            end
        end
    end
    %Create one aditional association in kegg:
    n         = n+1;
    kegg{n,1} = uni;
    kegg{n,2} = prot_name;
    kegg{n,3} = gene_name;
    kegg{n,4} = EC_names;
    kegg{n,5} = MW;
    kegg{n,6} = pathway;
    kegg{n,7} = sequence;
    
    %Update complex database:
    [kegg_complex,m] = updateComplexDB(kegg_complex,m,prot_name,uni,MW);
    cd ../../Databases/KEGG
    disp(['Updating KEGG database: Ready with gene ' gene_name])
end
kegg(n+1:end,:)         = [];
kegg_complex(m+1:end,:) = [];

%Save all databases as .mat files:
cd ..
save('ProtDatabase.mat','kegg','kegg_complex','swissprot','swissprot_complex');
cd ../Matlab_Module/get_enzyme_data

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [complex,m] = updateComplexDB(complex,m,prot_name,uniprot,MW)
%Adds protein to the complex database, as long as a "subunit" appears in
%the protein name.

subunit_pos = strfind(prot_name,'subunit');
if ~isempty(subunit_pos)
    complex_name = prot_name(1:subunit_pos-2);
    match        = false;
    for j = 1:m
        if strcmp(complex_name,complex{j,1})
            %If protein is part of a complex, update complex database:
            match        = true;
            complex{j,2} = [complex{j,2} '-' uniprot];
            complex{j,3} = complex{j,3} + MW;
        end
    end
    %If no match, create new complex:
    if ~match
        m            = m+1;
        complex{m,1} = complex_name;
        complex{m,2} = uniprot;
        complex{m,3} = MW;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%