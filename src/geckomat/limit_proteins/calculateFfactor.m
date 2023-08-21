function f = calculateFfactor(model, protData, enzymes, modelAdapter)
% calculateFfactor
%   Computes the f factor, as a proxy to the mass fraction of proteins
%   accounted for in an ecModel out of the total protein content in cells.
%
% Input:
%   model        an ecModel in GECKO 3 format (with ecModel.ec structure)
%   protData     structure with proteome data, from loadProtData (Optional,
%                by default it instead attempts to load data/paxDB.tsv)
%   enzymes      list of enzymes (Optional, default model.ec.enzymes)
%   modelAdapter a loaded model adapter (Optional, will otherwise use the
%                default model adapter).
%
% Output:
%   f            f-factor
%
% Usage:
%   f = calculateFfactor(model, protData, enzymes, modelAdapter)

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 3 || isempty(enzymes)
    enzymes = model.ec.enzymes;
end

% Gather proteome data in protData structure
if nargin < 2 || isempty(protData)
    if exist(fullfile(params.path,'data','paxDB.tsv'),'file')
        protData = fullfile(params.path,'data','paxDB.tsv');
    else
        printOrange('WARNING: No proteomics data is provided or can be found. Default f value of 0.5 is returned.\n')
        f = 0.5;
    end
end

% Gather Uniprot database for finding MW
uniprotDB = loadDatabases('uniprot', modelAdapter);
uniprotDB = uniprotDB.uniprot;

if ischar(protData) && endsWith(protData,'paxDB.tsv')
    fID         = fopen(fullfile(protData),'r');
    fileContent = textscan(fID,'%s','delimiter','\n');
    headerLines = sum(startsWith(fileContent{1},'#'));
    fclose(fID);

    %Read data file, excluding headerlines
    fID         = fopen(fullfile(protData),'r');
    fileContent = textscan(fID,'%s %s %f','delimiter','\t','HeaderLines',headerLines);
    genes       = fileContent{2};
    %Remove internal geneIDs modifiers
    genes       = regexprep(genes,'^\d+\.','');
    level       = fileContent{3};
    fclose(fID);
    [a,b]       = ismember(genes,uniprotDB.genes);
    uniprot     = uniprotDB.ID(b(a));
    level(~a)   = [];
    clear protData
    protData.uniprot = uniprot;
    protData.level   = level;
end

% Get MW and abundance (unit does not matter, f is fraction)
[~,idx] = ismember(protData.uniprot,uniprotDB.ID);
protData.MW = uniprotDB.MW(idx);
protData.abundance = protData.level .* protData.MW;

totalProt = sum(protData.abundance);

% Get enzymes in model
enzymesInModel = ismember(protData.uniprot,enzymes);
totalEnz = sum(protData.abundance(enzymesInModel));

f = totalEnz/totalProt;
end
