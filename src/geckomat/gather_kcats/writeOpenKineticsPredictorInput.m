function writeOpenKineticsPredictorInput(model, ecRxns, modelAdapter, onlyWithSmiles, overwrite)
% writeOpenKineticsPredictorInput
%   Prepares input for OpenKineticsPredictor: protein sequences and single-
%   substrate SMILES, which allows the use of the following predictors:
%   DLKcat; EITLEM; UniKP; KinForm-H; KinForm-L; CataPro; CatPred.
% 
%   The function writes OKP.csv, that can be used as input file at
%   https://predictor.openkinetics.org/
%
%   The OpenKineticsPrecitor output file can be used to generate a kcatList
%   through readOpenKineticsPredictorOutput.
%
% Input:
%   model           ecModel in GECKO 3 format (with ecModel.ec structure)
%   ecRxns          logical vector indicating which reactions to include
%                   (Optional, default: all reactions)
%   modelAdapter    loaded model adapter (Optional, uses default if empty)
%   onlyWithSmiles  logical, only include entries with SMILES (Optional, default: true)
%   overwrite       logical, overwrite existing file (Optional, default: false)
%
% Usage:
%   writeOpenKineticsPredictorInput(model, ecRxns)

[geckoPath, ~] = findGECKOroot();

%% Parse inputs
if nargin<2 || isempty(ecRxns)
    ecRxns = true(numel(model.ec.rxns),1);
elseif ~islogical(ecRxns) || numel(ecRxns) ~= numel(model.ec.rxns)
    error('ecRxns should be a logical vector with length equal to model.ec.rxns')
end

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either provide a modelAdapter or set the default in ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin<4 || isempty(onlyWithSmiles)
    onlyWithSmiles = true;
end

if nargin<5 || isempty(overwrite)
    overwrite = false;
end

%% Output files
filename = fullfile(params.path,'data','OKP.csv');
if ~overwrite && exist(filename,'file')
    error('%s already exists. Set overwrite=true to replace it.', filename)
end

%% Map ec.rxns to model.rxns
if model.ec.geckoLight
   origRxns = extractAfter(model.ec.rxns, 4);
else
   origRxns = model.ec.rxns;
end

ecRxnsIdx = find(ecRxns);
[found, origRxnIdxs] = ismember(origRxns(ecRxnsIdx), model.rxns);
if ~all(found)
    error('Not all reactions in model.ec.rxns found in model.rxns')
end

%% Build reduced stoichiometric matrix (ignore certain metabolites)
% Load metabolites to ignore (reuses DLKcat files)
if exist(fullfile(params.path,'data','DLKcatIgnoreMets.tsv'),'file')
    ignoreFile = fullfile(params.path,'data','DLKcatIgnoreMets.tsv');
else
    ignoreFile = fullfile(geckoPath,'databases','DLKcatIgnoreMets.tsv');
end
fID = fopen(ignoreFile);
fileData = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[ignoreMets, ignoreSmiles] = deal(fileData{[1,2]});

% Normalize metabolite names (remove special characters, lowercase)
metsNorm        = lower(regexprep(model.metNames,'[^0-9a-zA-Z]+',''));
ignoreMetsNorm  = lower(regexprep(ignoreMets,'[^0-9a-zA-Z]+',''));
ignoreSmiles(cellfun(@isempty,ignoreSmiles)) = [];

% Mark metabolites to ignore
ignoreMetsIdx = ismember(metsNorm, ignoreMetsNorm);
if isfield(model,'metSmiles')
    ignoreMetsIdx = ignoreMetsIdx | ismember(model.metSmiles, ignoreSmiles);
end
ignoreMetsIdx = ignoreMetsIdx | startsWith(model.mets,'prot_');

reducedS = model.S;
reducedS(ignoreMetsIdx,:) = 0;

%% Remove currency metabolite pairs (reuses DLKcat files)
if exist(fullfile(params.path,'data','DLKcatCurrencyMets.tsv'),'file')
    currencyFile = fullfile(params.path,'data','DLKcatCurrencyMets.tsv');
else
    currencyFile = fullfile(geckoPath,'databases','DLKcatCurrencyMets.tsv');
end
fID = fopen(currencyFile);
fileData = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[currencyMets(:,1), currencyMets(:,2)] = deal(fileData{[1,2]});
currencyMets = lower(regexprep(currencyMets,'[^0-9a-zA-Z]+',''));
for i=1:size(currencyMets,1)
    subs = strcmp(currencyMets(i,1),metsNorm);
    prod = strcmp(currencyMets(i,2),metsNorm);
    [~,subsRxns]=find(reducedS(subs,:));
    [~,prodRxns]=find(reducedS(prod,:));
    pairRxns = intersect(subsRxns,prodRxns);
    tempRedS=reducedS;
    tempRedS([find(subs);find(prod)],pairRxns) = 0;
    % Do not remove currency mets if no substrate remains
    rxnsWithRemainingSubstrates = any(tempRedS(:,pairRxns) < 0,1);
    reducedS([find(subs);find(prod)],intersect(pairRxns,pairRxns(rxnsWithRemainingSubstrates))) = 0;
end


%% Extract substrates for selected reactions only
clearedS = reducedS(:, origRxnIdxs);

% Find all substrates (negative stoichiometry)
[substrateIdxs, reactionIdxs] = find(clearedS < 0);

% Find proteins catalyzing each reaction
[proteinIdxs, rxnEnzIdxs] = find(model.ec.rxnEnzMat(reactionIdxs,:)');

%% Build output: sequence and SMILES
sequences = model.ec.sequence(proteinIdxs);

if isfield(model,'metSmiles')
    smiles = model.metSmiles(substrateIdxs(rxnEnzIdxs));
else
    smiles = cell(size(sequences));
end

%% Filter empty entries
if onlyWithSmiles
    valid = ~cellfun(@isempty, sequences) & ~cellfun(@isempty, smiles);
else
    valid = ~cellfun(@isempty, sequences);
    smiles(cellfun(@isempty, smiles)) = {'None'};
end

sequences = sequences(valid);
smiles = smiles(valid);

%% Write output file
% Sequence and SMILES
outTable = [sequences(:)'; smiles(:)'];

fID = fopen(filename,'w');
fprintf(fID,'Protein Sequence,Substrate\n');
fprintf(fID,'%s,%s\n', outTable{:});
fclose(fID);
fprintf('OpenKineticsPredictor input stored at %s\n', filename);
end
