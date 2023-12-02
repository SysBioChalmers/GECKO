function databases = loadDatabases(selectDatabase,modelAdapter)
% loadDatabases
%   Loads (and downloads if necessary) the organism-specific KEGG and
%   UniProt databases that are required to extract protein information. The
%   uniprot.ID and kegg.ID are taken from the ModelAdapter.
%
% Input:
%   selectDatabase  which databases should be loaded, either 'uniprot',
%                   'kegg' or 'both' (optional, default 'both')
%   modelAdapter    Model adapter. Optional, default will use the default 
%                   model adapter (send in [] for default).
%
% Output:
%   databases       contains .uniprot and .kegg structures, dependent on
%                   which databases were selected.
%
% Usage:
%   databases = loadDatabases(selectDatabase,modelAdapter)

if nargin<1
    selectDatabase = 'both';
end

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

params      = modelAdapter.getParameters();
kegg.ID      = params.kegg.ID;
uniprot.ID   = params.uniprot.ID;
filePath    = fullfile(params.path,'data');
uniprot.geneIDfield = params.uniprot.geneIDfield;
uniprot.type = params.uniprot.type;
if strcmp(uniprot.type,'taxonomy')
    uniprot.type = 'taxonomy_id';
end
kegg.geneID = params.kegg.geneID;
if params.uniprot.reviewed
    uniprotRev = 'reviewed:true+AND+';
else
    uniprotRev = '';
end

warning('off', 'MATLAB:MKDIR:DirectoryExists');

databases.uniprot = [];
databases.kegg = [];

%% Uniprot
if any(strcmp(selectDatabase,{'uniprot','both'}))
    uniprotPath = fullfile(filePath,'uniprot.tsv');
    if ~exist(uniprotPath,'file')
        if isempty(uniprot.ID)
            printOrange('WARNING: No uniprot.ID is specified, unable to download UniProt DB.\n')
        end
        disp(['Downloading Uniprot data for ' uniprot.type ' ' uniprot.ID '. This can take a few minutes.'])
        url = ['https://rest.uniprot.org/uniprotkb/stream?query=' uniprotRev ...
               uniprot.type ':' num2str(uniprot.ID) '&fields=accession%2C' uniprot.geneIDfield ...
            '%2Cec%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
        try
            urlwrite(url,uniprotPath,'Timeout',30);
            fprintf('Model-specific KEGG database stored at %s\n',uniprotPath);
        catch
            error(['Download failed, check your internet connection and try again, or manually download: ' url])
        end
    end
    if exist(uniprotPath,'file')
        fid         = fopen(uniprotPath,'r');
        fileContent = textscan(fid,'%q %q %q %q %q','Delimiter','\t','HeaderLines',1);
        fclose(fid);
        databases.uniprot.ID      = fileContent{1};
        databases.uniprot.genes   = fileContent{2};
        databases.uniprot.eccodes = fileContent{3};
        databases.uniprot.MW      = str2double(fileContent{4});
        databases.uniprot.seq     = fileContent{5};
    else
        databases.uniprot = [];
    end
    if ~isempty(databases.uniprot)
        [uniqueIDs,uniqueIdx] = unique(databases.uniprot.ID,'stable');
        if numel(uniqueIDs) < numel(databases.uniprot.ID)
            duplID = setdiff(1:numel(databases.uniprot.ID),uniqueIdx);
            dispEM(['Duplicate entries are found for the following proteins. '...
                    'Manually curate the ''uniprot.tsv'' file, or adjust the uniprot parameters '...
                    'in the model adapter:'],true,databases.uniprot.ID(duplID));
        end
    end
end

%% KEGG
if any(strcmp(selectDatabase,{'kegg','both'}))
    keggPath = fullfile(filePath,'kegg.tsv');
    if ~exist(keggPath,'file')
        if isempty(kegg.ID)
            printOrange('WARNING: No kegg.ID is specified, unable to download KEGG DB.\n')
        else
            downloadKEGG(kegg.ID,keggPath,kegg.geneID);
        end
    end
    if exist(keggPath,'file')
        fid         = fopen(keggPath,'r');
        fileContent = textscan(fid,'%q %q %q %q %q %q %q','Delimiter',',','HeaderLines',0);
        fclose(fid);
        databases.kegg.uniprot    = fileContent{1};
        databases.kegg.genes      = fileContent{2};
        databases.kegg.keggGene   = fileContent{3};
        databases.kegg.eccodes    = fileContent{4};
        databases.kegg.MW         = str2double(fileContent{5});
        databases.kegg.pathway    = fileContent{6};
        databases.kegg.seq        = fileContent{7};
    else
        databases.kegg = [];
    end
end
end

function downloadKEGG(keggID, filePath, keggGeneID)
%% Download gene information
progressbar(['Downloading KEGG data for organism code ' keggID])
webOptions = weboptions('Timeout',30);
try
    gene_list = webread(['http://rest.kegg.jp/list/' keggID],webOptions);
catch ME
    switch ME.identifier
        case 'MATLAB:webservices:HTTP400StatusCodeError'
            error(['Unable to download data form KEGG with a potentially invalid ID: ' keggID ])
    end
end
gene_list = regexpi(gene_list, '[^\n]+','match')';
gene_id   = regexpi(gene_list,['(?<=' keggID ':)\S+'],'match');

% Retrieve information for every gene in the list, 10 genes per query
genesPerQuery = 10;
queries = ceil(numel(gene_id)/genesPerQuery);
keggData  = cell(numel(gene_id),1);
for i = 1:queries
    % Download batches of genes
    firstIdx = i*genesPerQuery-(genesPerQuery-1);
    lastIdx  = i*genesPerQuery;
    if lastIdx > numel(gene_id) % Last query has probably less genes
        lastIdx = numel(gene_id);
    end
    url      = ['http://rest.kegg.jp/get/' keggID ':' strjoin([gene_id{firstIdx:lastIdx}],['+' keggID ':'])];

    retry = true;
    while retry
        try
            retry = false;
            out   = webread(url,webOptions);
        catch
            retry = true;
        end
    end
    outSplit = strsplit(out,['///' 10]); %10 is new line character
    if numel(outSplit) < lastIdx-firstIdx+2
        error('KEGG returns less genes per query') %Reduce genesPerQuery
    end
    keggData(firstIdx:lastIdx) = outSplit(1:end-1);
    progressbar(i/queries)
end

%% Parsing of info to keggDB format
sequence  = regexprep(keggData,'.*AASEQ\s+\d+\s+([A-Z\s])+?\s+NTSEQ.*','$1');
%No AASEQ -> no protein -> not of interest
noProt    = startsWith(sequence,'ENTRY ');
uni       = regexprep(keggData,'.*UniProt: (\S+?)\s.*','$1');
noUni     = startsWith(uni,'ENTRY ');
uni(noProt | noUni)       = [];
keggData(noProt | noUni) = [];
sequence(noProt | noUni)  = [];
sequence  = regexprep(sequence,'\s*','');
keggGene  = regexprep(keggData,'ENTRY\s+(\S+?)\s.+','$1');

switch keggGeneID
    case {'kegg',''}
        gene_name = keggGene;
    otherwise
        % In case there are special characters:
        keggGeneIDT = regexptranslate('escape',keggGeneID);
        gene_name = regexprep(keggData,['.+' keggGeneIDT ': (\S+?)\n.+'],'$1');
        noID = ~contains(keggData,keggGeneIDT);
        if all(noID)
            error(['None of the KEGG entries are annotated with the gene identifier ' keggGeneID])
        else
            gene_name(noID)= [];
            keggData(noID) = [];
            keggGene(noID) = [];
            sequence(noID) = [];
            uni(noID)      = [];
        end
end

EC_names  = regexprep(keggData,'.*ORTHOLOGY.*\[EC:(.*?)\].*','$1');
EC_names(startsWith(EC_names,'ENTRY ')) = {''};

MW = cell(numel(sequence),1);
for i=1:numel(sequence)
    if ~isempty(sequence{i})
        MW{i} = num2str(round(calculateMW(sequence{i})));
    end
end

pathway   = regexprep(keggData,'.*PATHWAY\s+(.*?)(BRITE|MODULE).*','$1');
pathway(startsWith(pathway,'ENTRY ')) = {''};
pathway   = strrep(pathway,[keggID '01100  Metabolic pathways'],'');
pathway   = regexprep(pathway,'\n','');
pathway   = regexprep(pathway,'           ','');

out = [uni, gene_name, keggGene, EC_names, MW, pathway, sequence];
out = cell2table(out);

writetable(out, filePath, 'FileType', 'text', 'WriteVariableNames',false);
fprintf('Model-specific KEGG database stored at %s\n',filePath);
end
