function f = calculateFfactor(model, varargin)
% calculateFfactor  Compute the f factor for an ecModel.
%
% Computes the f factor, as a proxy to the mass fraction of proteins
% accounted for in an ecModel out of the total protein content in cells.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% protData : struct
%     structure with proteome data, from loadProtData (by default it instead
%     attempts to load data/paxDB.tsv).
% enzymes : cell
%     list of enzymes (default model.ec.enzymes).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% f : double
%     f-factor.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     f = calculateFfactor(model);
%     f = calculateFfactor(model, 'enzymes', enzymes);
%
% See also
% --------
% loadProtData, getProtFromProteomics

p = parseGECKOargs(varargin, { ...
    'protData',     []; ...
    'enzymes',      []; ...
    'modelAdapter', []});
protData     = p.protData;
enzymes      = p.enzymes;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(enzymes)
    enzymes = model.ec.enzymes;
end

% Gather proteome data in protData structure
if isempty(protData)
    if exist(fullfile(params.path,'data','paxDB.tsv'),'file')
        protData = fullfile(params.path,'data','paxDB.tsv');
    else
        printOrange('WARNING: No proteomics data is provided or can be found. Default f value of 0.5 is returned.\n');
        f = 0.5;
        return
    end
end

if ischar(protData) && endsWith(protData,'paxDB.tsv')
    % Gather Uniprot database for finding MW. Only needed here, to map paxDB.tsv's gene
    % identifiers onto UniProt IDs and molecular weights --- a protData struct supplied by
    % the caller already carries UniProt IDs and needs no lookup, so this stays inside the
    % branch that actually uses it rather than running (and potentially downloading or
    % reading UniProt data) on every call.
    uniprotDB = loadDatabases('uniprot', modelAdapter);
    uniprotDB = uniprotDB.uniprot;

    fID         = fopen(fullfile(protData),'r');
    fileContent = textscan(fID,'%s','delimiter','\n');
    headerLines = find(startsWith(fileContent{1},'#'),1,'last');
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
    protData.uniprotIDs = uniprot;
    protData.level   = level;
    % Get MW and abundance (unit does not matter, f is fraction)
    [~,idx] = ismember(protData.uniprotIDs,uniprotDB.ID);
    protData.MW = uniprotDB.MW(idx);
    protData.abundances = protData.level .* protData.MW;
end

avgAbundances = mean(protData.abundances,2);
totalProt = sum(avgAbundances,'omitnan');

% Get enzymes in model
enzymesInModel = ismember(protData.uniprotIDs,enzymes);
totalEnz = sum(avgAbundances(enzymesInModel),'omitnan');

f = totalEnz/totalProt;
end
