function kcatList = readOpenKineticsPredictorOutput(model, outFile, modelAdapter)
% readOpenKineticsPredictorOutput
%   Reads the OpenKineticsPredictor output file (job-*-output.csv with
%   columns: kcat (1/s), Source kcat, Extra Info kcat, Protein Sequence,
%   Substrate) and constructs a kcatList structure, that can be used by
%   selectKcatValue() to populate the ecModel with kcat values.
%
%   The per-entry `Source kcat` column (e.g. 'Prediction from CatPred',
%   'BRENDA', 'Sabio-RK', 'UniProt') is preserved in kcatList.kcatSource
%   so that selectKcatValue() records the actual provenance of each value.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   outFile         name and path of the OpenKineticsPredictor output CSV.
%                   (Optional; if omitted, a file selection dialog is
%                   opened, starting in the data/ folder specified by the
%                   modelAdapter)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
% Output:
%   kcatList    structure array with list of kcat values from
%               OpenKineticsPredictor, with separate entries for each kcat
%               value
%               source      'OpenKineticsPredictor'
%               rxns        reaction identifiers, matching model.ec.rxns
%               genes       gene identifiers, matching model.ec.genes
%               substrates  substrate names, matching model.metNames
%               kcats       kcat values in /sec
%               kcatSource  per-entry source as reported by OKP, prefixed
%                           with 'OKP-' (e.g. 'OKP-CatPred', 'OKP-BRENDA',
%                           'OKP-Sabio-RK', 'OKP-UniProt')
%
% Usage:
%   kcatList = readOpenKineticsPredictorOutput(model, outFile, modelAdapter)

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin < 2 || isempty(outFile)
    startDir = fullfile(params.path,'data');
    if ~exist(startDir,'dir')
        startDir = params.path;
    end
    [fName, fPath] = uigetfile({'*.csv','OpenKineticsPredictor output (*.csv)'; ...
                                '*.*','All files (*.*)'}, ...
                               'Select OpenKineticsPredictor output file', ...
                               fullfile(startDir, filesep));
    if isequal(fName,0)
        error('No file selected.')
    end
    outFile = fullfile(fPath, fName);
end

%% Read OKP output file using readcell to correctly handle quoted fields
% that may contain commas in the 'Extra Info kcat' column.
if ~exist(outFile,'file')
    error('Cannot find file: %s', outFile)
end
raw = readcell(outFile, 'Delimiter', ',', 'TextType', 'char');

if size(raw,2) < 5
    error('OKP output file does not have the expected 5 columns (kcat, Source kcat, Extra Info kcat, Protein Sequence, Substrate): %s', outFile)
end

% Drop header row
raw = raw(2:end, :);

kcatsCol      = raw(:,1);
sourceKcatCol = raw(:,2);
sequences     = raw(:,4);
smiles        = raw(:,5);

% Normalize types: readcell returns doubles for numeric cells, char for
% strings, and `missing` for empty fields.
kcats = nan(size(kcatsCol));
for i = 1:numel(kcatsCol)
    v = kcatsCol{i};
    if isnumeric(v)
        kcats(i) = v;
    elseif ischar(v) || isstring(v)
        kcats(i) = str2double(v);
    end
end

sourceKcatCol = normalizeCellStrings(sourceKcatCol);
sequences     = normalizeCellStrings(sequences);
smiles        = normalizeCellStrings(smiles);

%% Filter out entries with no numeric kcat value (NaN). Zero values are
% retained here and filtered downstream by selectKcatValue.
valid = ~isnan(kcats);
kcats         = kcats(valid);
sequences     = sequences(valid);
smiles        = smiles(valid);
sourceKcatCol = sourceKcatCol(valid);

if isempty(kcats)
    error('OKP file does not contain any numeric kcat values. Run the job at https://predictor.openkinetics.org/ first.')
end

if ~isfield(model,'metSmiles')
    error('model.metSmiles is required to map SMILES back to metabolites.')
end

%% Build ec.rxns -> model.rxns index mapping (inverse of the mapping done
%  in writeOpenKineticsPredictorInput)
if model.ec.geckoLight
    origRxns = extractAfter(model.ec.rxns, 4);
else
    origRxns = model.ec.rxns;
end
[foundEcRxn, ecRxnToModelRxn] = ismember(origRxns, model.rxns);
if ~all(foundEcRxn)
    error('Not all entries in model.ec.rxns could be matched to model.rxns.')
end

%% Resolve each row to matching ec.rxns indices
rxnsOut       = cell(0,1);
genesOut      = cell(0,1);
substratesOut = cell(0,1);
kcatsOut      = zeros(0,1);
sourceOut     = cell(0,1);
unmatched     = 0;

for i = 1:numel(sequences)
    % Map protein sequence -> ec.genes indices (columns of rxnEnzMat)
    proteinIdxs = find(strcmp(model.ec.sequence, sequences{i}));
    % Map SMILES -> metabolite indices (rows of model.S)
    metIdxs = find(strcmp(model.metSmiles, smiles{i}));

    if isempty(proteinIdxs) || isempty(metIdxs)
        unmatched = unmatched + 1;
        continue
    end

    % model.rxns columns where any of these mets are substrates (S < 0)
    modelRxnsUsingSub = any(model.S(metIdxs, :) < 0, 1);

    % Lift to ec.rxns: does the underlying model.rxn consume this metabolite?
    ecRxnUsesSub = modelRxnsUsingSub(ecRxnToModelRxn).';

    % ec.rxns catalyzed by any of these protein(s)
    ecRxnCatalyzed = any(model.ec.rxnEnzMat(:, proteinIdxs) ~= 0, 2);

    candidateEcRxns = find(ecRxnUsesSub & ecRxnCatalyzed);

    if isempty(candidateEcRxns)
        unmatched = unmatched + 1;
        continue
    end

    for k = 1:numel(candidateEcRxns)
        ecIdx = candidateEcRxns(k);

        % Pick the gene whose protein both matches this sequence and
        % catalyzes this ec.rxn
        catalysts       = find(model.ec.rxnEnzMat(ecIdx, :) ~= 0);
        matchingProtein = intersect(proteinIdxs, catalysts, 'stable');
        if ~isempty(matchingProtein)
            geneName = model.ec.genes{matchingProtein(1)};
        else
            geneName = '';
        end

        % Pick a metabolite name whose SMILES matches and that is a
        % substrate of the underlying model.rxn
        origRxnIdx       = ecRxnToModelRxn(ecIdx);
        metIsSub         = model.S(metIdxs, origRxnIdx) < 0;
        substrateMetIdxs = metIdxs(metIsSub);
        if ~isempty(substrateMetIdxs)
            subName = model.metNames{substrateMetIdxs(1)};
        else
            subName = '';
        end

        rxnsOut{end+1,1}       = model.ec.rxns{ecIdx};
        genesOut{end+1,1}      = geneName;
        substratesOut{end+1,1} = subName;
        kcatsOut(end+1,1)      = kcats(i);
        sourceOut{end+1,1}     = sourceKcatCol{i};
    end
end

if unmatched > 0
    warning('%d of %d OKP entries could not be matched to an ec.rxn (protein sequence and/or SMILES not found in model, or no enzyme+substrate pair matches).', ...
        unmatched, numel(sequences));
end

if isempty(rxnsOut)
    error('No OKP entries could be mapped to model.ec.rxns.')
end

%% Build per-entry source labels
% Strip leading 'Prediction from ' to keep labels compact, then prefix all
% with 'OKP-' so the overall pipeline stage is still identifiable in
% model.ec.source.
kcatSource = regexprep(sourceOut, '^Prediction from\s+', '');
kcatSource = strcat('OKP-', kcatSource);

%% Make kcatList structure
kcatList.source     = 'OpenKineticsPredictor';
kcatList.rxns       = rxnsOut;
kcatList.genes      = genesOut;
kcatList.substrates = substratesOut;
kcatList.kcats      = kcatsOut;
kcatList.kcatSource = kcatSource;
end

function c = normalizeCellStrings(c)
% Convert each cell entry to a char row vector; missing/empty -> ''.
for i = 1:numel(c)
    v = c{i};
    if ismissing(v)
        c{i} = '';
    elseif isstring(v)
        c{i} = char(v);
    elseif isnumeric(v)
        if isnan(v)
            c{i} = '';
        else
            c{i} = num2str(v);
        end
    elseif ~ischar(v)
        c{i} = char(string(v));
    end
end
end
