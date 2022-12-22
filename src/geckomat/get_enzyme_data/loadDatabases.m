function databases = loadDatabases(selectDatabase)
% loadDatabases
%   Load or downloads (if necessary) the organism-specific KEGG and UniProt
%   databases that are required to extract protein information.
%
% Input:
%   selectDatabase  which databases should be loaded, either 'uniprot',
%                   'kegg' or 'both'
%
% Output:
%   databases       contains .uniprot and .kegg structures, dependent on
%                   which databases were selected.
%
% Usage: databases = loadDatabases(selectDatabase)

global GECKOModelAdapter
params=checkGECKOModelAdapter(GECKOModelAdapter);
keggID=param.keggID;
taxonID=param.keggID;

if nargin<1
    selectDatabase = 'both';
end

weboptions('Timeout',10);
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Uniprot
if any(strcmp(selectDatabase,{'uniprot','both'}))
    filePath = fullfile(params.path,'data','uniprot.tsv');
    if ~exist(filePath,'file')
        if isempty(taxonID)
            error('No taxonID is specified, unable to download UniProt DB')
        end
        disp(['Downloading Uniprot data for taxonomic ID ' taxonID '. This can take a few minutes.'])
        url = ['https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:' num2str(taxonID) ...
            '&fields=accession%2C' param.uniprotGeneIdField ...
            '%2Cgene_primary%2Cec%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
        urlwrite(url,filePath);
    end

    fid         = fopen(filePath,'r');
    fileContent = textscan(fid,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
    fclose(fid);
    databases.uniprot.ID      = fileContent{1};
    databases.uniprot.genes   = fileContent{2};
    databases.uniprot.eccodes = fileContent{4};
    databases.uniprot.MW      = str2double(fileContent{5});
    databases.uniprot.seq     = fileContent{6};
end

%% KEGG
if any(strcmp(selectDatabase,{'kegg','both'}))
    filePath = fullfile(params.path,'data','kegg.tsv');
    if ~exist(filePath,'file')
        if isempty(keggID)
            error('No keggID is specified, unable to download KEGG DB')
        end
        downloadKEGG(keggID,filePath,verbose);
    end

    fid         = fopen(filePath,'r');
    fileContent = textscan(fid,'%s %s %s %s %s %s','Delimiter',',','HeaderLines',0);
    fclose(fid);

    databases.kegg.uniprot    = fileContent{1};
    databases.kegg.genes      = fileContent{2};
    databases.kegg.eccodes    = fileContent{3};
    databases.kegg.MW         = str2double(fileContent{4});
    databases.kegg.pathway    = fileContent{5};
    databases.kegg.seq        = fileContent{6};
end
end

function downloadKEGG(keggID, filePath, verbose)
%% Download gene information
if verbose
    fprintf('Downloading KEGG data for organism code %s...   0%% complete',keggID);
end
gene_list = webread(['http://rest.kegg.jp/list/' keggID]);
gene_list = regexpi(gene_list, '[^\n]+','match')';
gene_id   = regexpi(gene_list,['(?<=' keggID ':)\S+'],'match');

% Retrieve information for every gene in the list, 10 genes per query
genesPerQuery = 10;
queries = ceil(numel(gene_id)/genesPerQuery);
keggGenes  = cell(numel(gene_id),1);
for i = 1:queries
    % Report progress
    if verbose
        progress=num2str(floor(100*i/queries));
        progress=pad(progress,3,'left');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
    end
    % Download batches of genes
    firstIdx = i*genesPerQuery-(genesPerQuery-1);
    lastIdx  = i*genesPerQuery;
    if lastIdx > numel(gene_id) % Last query has probably less genes
        lastIdx = numel(gene_id);
    end
    url      = ['http://rest.kegg.jp/get/' keggID ':' strjoin([gene_id{firstIdx:lastIdx}],['+' keggID ':'])];

    out      = webread(url);
    outSplit = strsplit(out,['///' 10]); %10 is new line character
    if numel(outSplit) < lastIdx-firstIdx+2
        error('KEGG returns less genes per query') %Reduce genesPerQuery
    end
    keggGenes(firstIdx:lastIdx) = outSplit(1:end-1);
end

%% Parsing of info to keggDB format
sequence  = regexprep(keggGenes,'.*AASEQ\s+\d+\s+([A-Z\s])+?\s+NTSEQ.*','$1');
%No AASEQ -> no protein -> not of interest
noProt    = startsWith(sequence,'ENTRY ');
keggGenes(noProt) = [];
sequence(noProt) = [];
sequence  = regexprep(sequence,'\s*','');

uni       = regexprep(keggGenes,'.*UniProt: (\S+?)\s.*','$1');
uni(startsWith(uni,'ENTRY ')) = {''};

gene_name = regexprep(keggGenes,'ENTRY\s+(\S+?)\s.+','$1');

% prot_name = regexprep(keggData,'.*NAME\s+(\S.+?)\n.+','$1');
% prot_name = regexprep(prot_name,'\(RefSeq\) ','');

EC_names  = regexprep(keggGenes,'.*ORTHOLOGY.*\[EC:(.*?)\].*','$1');
EC_names(startsWith(EC_names,'ENTRY ')) = {''};

sequence  = regexprep(keggGenes,'.*AASEQ\s+\d+\s+([A-Z\s])+?\s+NTSEQ.*','$1');
sequence(startsWith(sequence,'ENTRY ')) = {''};
sequence  = regexprep(sequence,'\s*','');

MW = cell(numel(sequence),1);
for i=1:numel(sequence)
    if ~isempty(sequence{i})
        MW{i} = num2str(round(calculateMW(sequence{i})));
    end
end

pathway   = regexprep(keggGenes,'.*PATHWAY\s+(.*?)(BRITE|MODULE).*','$1');
pathway(startsWith(pathway,'ENTRY ')) = {''};
pathway   = strrep(pathway,[keggID '01100  Metabolic pathways'],'');
pathway   = regexprep(pathway,'\n','');
pathway   = regexprep(pathway,'           ','');

out = [uni, gene_name, EC_names, MW, pathway, sequence];
out = cell2table(out);

writetable(out, filePath, 'FileType', 'text', 'WriteVariableNames',false);
if verbose
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b100%% complete\n');
end
end
