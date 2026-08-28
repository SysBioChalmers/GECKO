function databases = loadDatabases(varargin)
% loadDatabases  Load the organism-specific KEGG and UniProt databases.
%
% Loads (and downloads if necessary) the organism-specific KEGG and UniProt
% databases that are required to extract protein information. The uniprot.ID
% and kegg.ID are taken from the ModelAdapter.
%
% Name-Value Arguments
% --------------------
% selectDatabase : char
%     which databases should be loaded, either 'uniprot', 'kegg' or 'both'
%     (default 'both').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter;
%     send in [] for default).
%
% Returns
% -------
% databases : struct
%     contains .uniprot and .kegg structures, dependent on which databases
%     were selected.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     databases = loadDatabases('both', adapter);
%     databases = loadDatabases('selectDatabase', 'both');
%
% See also
% --------
% getECfromDatabase, loadBRENDAdata, downloadKEGG, downloadUniProt

p = parseGECKOargs(varargin, { ...
    'selectDatabase', 'both'; ...
    'modelAdapter',   []});
selectDatabase = p.selectDatabase;
modelAdapter   = p.modelAdapter;

if isempty(modelAdapter)
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
kegg.geneID = params.kegg.geneID;

warning('off', 'MATLAB:MKDIR:DirectoryExists');

databases.uniprot = [];
databases.kegg = [];

%% Uniprot
if any(strcmp(selectDatabase,{'uniprot','both'}))
    uniprotPath = fullfile(filePath,'uniprot.tsv');
    if ~exist(uniprotPath,'file')
        if isempty(uniprot.ID)
            printOrange('WARNING: No uniprot.ID is specified, unable to download UniProt DB.\n');
        end
        downloadUniProt(uniprot.ID, uniprotPath, uniprot.geneIDfield, uniprot.type, params.uniprot.reviewed);
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
            error('%s', ravenList(['Duplicate entries are found for the following proteins. '...
                    'Manually curate the ''uniprot.tsv'' file, or adjust the uniprot parameters '...
                    'in the model adapter:'],databases.uniprot.ID(duplID)));
        end
    end
end

%% KEGG
if any(strcmp(selectDatabase,{'kegg','both'}))
    keggPath = fullfile(filePath,'kegg.tsv');
    if ~exist(keggPath,'file')
        if isempty(kegg.ID)
            printOrange('WARNING: No kegg.ID is specified, unable to download KEGG DB.\n');
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
