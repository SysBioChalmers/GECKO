function [f,count] = calculateFfactor(model, protData, modelAdapter)
% calculateFfactor
%
% Computes the f factor, as a proxy to the mass fraction of proteins accounted 
% for in an ecModel out of the total protein content in cells. An
% integrated quantitative proteomics dataset from the https://pax-db.org/
% database (stored in this toolbox as: 'GECKO/databases/prot_abundance.txt')
% is used as a comparison basis.
% 
% Usage: [f,count] = measureAbundance(enzymes)
%

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

% Gather proteome data in protData structure
if nargin < 2 || isempty(protData)
    if exist(fullfile(params.path,'data','paxDB.tsv'),'file')
        protData = fullfile(params.path,'data','paxDB.tsv');
    elseif exist(fullfile(params.path,'data','proteomics.tsv'),'file')
        protData = fullfile(params.path,'data','proteomics.tsv');
    else
        disep('No proteomics data is provided or can be found. Default f value of 0.5 is returned.')
        f = 0.5;
    end
end

% Gather Uniprot database for finding MW
uniprotDB = loadDatabases('uniprot', modelAdapter);
uniprotDB = uniprotDB.uniprot;

if endsWith(protData,'paxDB.tsv')
    fID         = fopen(fullfile(protData),'r');
    fileContent = textscan(fID,'%s','delimiter','\n');
    headerLines = sum(startsWith(fileContent{1},'#'));
    fclose(fID);

    %Read data file, excluding headerlines
    fID         = fopen(fullfile(protData),'r');
    fileContent = textscan(fID,'%s %s %f','delimiter','\t','HeaderLines',headerLines);
    genes       = fileContent{2};
    %Remove internal geneIDs modifiers
    genes       = regexprep(genes,'(\d{4}).','');
    level       = fileContent{3};
    fclose(fID);
    [a,b]       = ismember(genes,uniprotDB.genes);
    uniprot     = uniprotDB.ID(b(a));
    level(~a)   = [];
    clear protData
    protData.uniprot = uniprot;
    protData.level   = level;
else
    [~, protData] = readProteomics(model, protData);
end

% Get MW and abundance (unit does not matter, f is fraction)
[~,idx] = ismember(protData.uniprot,uniprotDB.ID);
protData.MW = uniprotDB.MW(idx);
protData.abundance = protData.level .* protData.MW;

totalProt = sum(protData.abundance);

% Get enzymes in model
enzymes = ismember(protData.uniprot,model.ec.enzymes);
totalEnz = sum(protData.abundance(enzymes));

f = totalEnz/totalProt;
end
