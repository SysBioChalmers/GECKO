function databases = loadDatabases(selectDatabase,modelAdapter)
% loadDatabases
%   Loads (and downloads if necessary) the organism-specific KEGG and
%   UniProt databases that are required to extract protein information. The
%   uniprot.ID and kegg.ID are taken from the ModelAdapter.
%
%   Downloads, when triggered, are dispatched to the helpers
%   downloadKEGG.m and downloadUniProt.m.
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
filePath    = fullfile(params.path,'data');

warning('off', 'MATLAB:MKDIR:DirectoryExists');

databases.uniprot = [];
databases.kegg = [];

%% Uniprot
if any(strcmp(selectDatabase,{'uniprot','both'}))
    uniprotPath = fullfile(filePath,'uniprot.tsv');
    if ~exist(uniprotPath,'file')
        if isempty(params.uniprot.ID)
            printOrange('WARNING: No uniprot.ID is specified, unable to download UniProt DB.\n');
        else
            downloadUniProt(params.uniprot.ID, uniprotPath, ...
                params.uniprot.type, params.uniprot.geneIDfield, ...
                params.uniprot.reviewed);
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
        if isempty(params.kegg.ID)
            printOrange('WARNING: No kegg.ID is specified, unable to download KEGG DB.\n');
        else
            downloadKEGG(params.kegg.ID,keggPath,params.kegg.geneID);
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
