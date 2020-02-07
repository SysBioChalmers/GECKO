function [swissprot,kegg] = updateDatabases
% updateDatabases
%   Updates all databases for protein matching (KEGG and Swiss-Prot).
%
%    Note: Before using this script, one should manually download from 
%          http://www.uniprot.org/uniprot a tab delimited file for the
%          desired organism with the following format:
%          Entry - Protein names - Gene names - EC number - Sequence
%          OBS: filter with the Swiss-Prot option
% 
%   Usage: [swissprot,kegg] = updateDatabases
%
% Ivan Domenzain. Last edited: 2020-01-25
%

current = pwd;
cd ..
parameters = getModelParameters;
cd (current)
if isfield(parameters,'keggID')
    keggID = parameters.keggID;
    %Download KEGG data:
    mkdir ../../databases/KEGG
    downloadKEGGdata(keggID) 
    %Build KEGG table
    kegg = buildKEGGtable(keggID);
    %Remove KEGG files for compliance of repository:
    delete ../../databases/KEGG/*.txt
    rmdir ../../databases/KEGG
else
    kegg = cell(1,7);
end
%Build Swissprot table:
swissprot = buildSWISSPROTtable;
%Save both databases as .mat files:
save('../../databases/ProtDatabase.mat','swissprot','kegg');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function swissprot = buildSWISSPROTtable
%Build Swissprot table (uniprot code - protein name - gene names - EC number - sequence):
fileID_uni        = fopen('../../databases/uniprot.tab');
swissprot         = textscan(fileID_uni,'%s %s %s %s %s','delimiter','\t');
swissprot         = [swissprot{1} swissprot{2} swissprot{3} swissprot{4} swissprot{5}];
swissprot(1,:)    = [];
fclose(fileID_uni);
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
end
disp('Building Swiss-Prot database.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function downloadKEGGdata(keggID)
base      = 'http://rest.kegg.jp/';
operation = 'list/';
gene_list = urlread([base operation keggID]);
gene_list = regexpi(gene_list, '[^\n]+','match')'; 
gene_id   = regexpi(gene_list,['(?<=' keggID ':)\S+'],'match');
% Retrieve information for every gene in the list (with a maximum of 10,000
% to avoid bulk downloads)
operation = 'get/';
for i = 1:min([numel(gene_id),10000])
    try
        gene = urlread([base operation keggID ':' gene_id{i}{1}]);
        fid  = fopen(['../../databases/KEGG/' gene_id{i}{1} '.txt'],'w');
        fprintf(fid,'%s',gene);
        fclose(fid);
        disp(['Downloading KEGG data for ' gene_id{i}{1}])
    catch    
        disp(['Cannot find ' gene_id{i}{1} ' in KEGG']);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kegg = buildKEGGtable(keggID)
%Build KEGG table (uniprot code - protein name - systematic gene name - EC number - MW - pathway - sequence):
file_names      = dir('../../databases/KEGG/');
file_names(1:2) = [];
kegg            = cell(100000,7);
n               = 0;
for i = 1:length(file_names)
    file_name = file_names(i).name;
    %3rd column: systematic gene name
    gene_name = file_name(1:end-4);
    %Retrieve all data as a cell with all rows:
    fID  = fopen(['../../databases/KEGG/' file_name]);
    text = textscan(fID,'%s','delimiter','\t');
    fclose(fID);
    text = text{1};
    
    uni       = '';
    prot_name = '';
    EC_names  = '';
    sequence  = '';
    MW        = 0;
    pathway   = '';
    for j = 1:length(text)
        line = text{j};
        if length(line) > 10
            %1st column: uniprot code
            if strcmp(line(1:8),'UniProt:')
                uni = line(10:end);
                
            %2nd column: protein name
            elseif strcmp(line(1:10),'DEFINITION')
                if strcmp(line(13:20),'(RefSeq)')
                    prot_name = lower(line(22:end));
                else
                    prot_name = lower(line(13:end));
                    disp([gene_name ': no RefSeq'])
                    pause
                end
            
            %4th column: EC number
            elseif strcmp(line(1:9),'ORTHOLOGY')
                pos_EC = strfind(line,'[EC:');
                if ~isempty(pos_EC)
                    EC_names = line(pos_EC+4:end-1);
                end
                
            %5th column and 7th column: MW & sequence
            elseif strcmp(line(1:5),'AASEQ')
                end_seq  = false;
                for k = j+1:length(text)
                    if length(text{k}) > 10
                        if strcmp(text{k}(1:5),'NTSEQ')
                            end_seq = true;
                        end
                    end
                    if ~end_seq
                        sequence = [sequence text{k}];
                    end
                end
                MW = calculateMW(sequence);
                
            %6th column: pathway
            elseif strcmp(line(1:7),'PATHWAY')
                start    = strfind(line,keggID);
                pathway  = line(start(1):end);
                end_path = false;
                for k = j+1:length(text)
                    nospace = strrep(text{k},[keggID '01100  Metabolic pathways'],'');
                    nospace = strrep(nospace,' ','');
                    if length(nospace) > 10
                        if strcmp(nospace(1:3),keggID) && ~end_path
                            start    = strfind(text{k},keggID);
                            pathway  = [pathway ' ' text{k}(start(1):end)];
                        else
                            end_path = true;
                        end
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
    disp(['Updating KEGG database: Ready with gene ' gene_name])
end
kegg(n+1:end,:) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
