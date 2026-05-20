function [done, kcatList] = fetchOpenKineticsPredictor(model, useStored, jobId, wait, pollInterval, modelAdapter)
% fetchOpenKineticsPredictor
%   Checks the status of an OpenKineticsPredictor job and, when it has
%   finished, downloads the result to data/OKP_output.csv and parses it
%   into a kcatList (consumable by selectKcatValue). Replaces the manual
%   download + readOpenKineticsPredictorOutput step.
%
%   The API key (only needed when useStored=false) is resolved from the
%   OKP_API_KEY environment variable or data/okpApiKey.txt; see
%   submitOpenKineticsPredictor for how to obtain and store a key.
%
% Input:
%   model           ecModel in GECKO 3 format (with ecModel.ec structure;
%                   model.metSmiles is required to map results back)
%   useStored       logical (Optional, default false). If true, skip the
%                   API entirely and parse the already-downloaded
%                   data/OKP_output.csv. If false, query the API.
%   jobId           OKP job id (Optional). If empty, read from
%                   data/OKP_job.txt (written by submitOpenKineticsPredictor).
%   wait            logical (Optional, default false). If true, poll until
%                   the job finishes; if false, report once and return.
%   pollInterval    seconds between polls when wait=true (Optional,
%                   default 30).
%   modelAdapter    loaded model adapter (Optional, uses default if empty)
%
% Output:
%   done            logical, true if a result was obtained (downloaded or
%                   loaded from disk)
%   kcatList        structure with kcat values (empty when not done):
%                   rxns; genes; substrates; kcats (1/s); kcatSource
%                   (per-entry actual provenance from the 'Source kcat'
%                   column: the prediction method such as 'CataPro', or a
%                   database such as 'BRENDA' / 'Sabio-RK' / 'UniProt');
%                   source (the most frequent of those, as a scalar label)
%
% Usage:
%   [done, kcatList] = fetchOpenKineticsPredictor(model);                  % check once
%   [done, kcatList] = fetchOpenKineticsPredictor(model, false, [], true); % poll
%   [done, kcatList] = fetchOpenKineticsPredictor(model, true);            % parse stored file

if nargin < 2 || isempty(useStored); useStored = false; end
if nargin < 6 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params  = modelAdapter.params;
if nargin < 4 || isempty(wait); wait = false; end
if nargin < 5 || isempty(pollInterval); pollInterval = 30; end

outFile = fullfile(params.path,'data','OKP_output.csv');

%% Stored-result mode: parse the saved file, no API call
if useStored
    if ~exist(outFile,'file')
        error('useStored=true but no stored result found at %s. Run fetchOpenKineticsPredictor with useStored=false first.', outFile);
    end
    kcatList = parseOkpOutput(model, outFile);
    done = true;
    return
end

%% API mode
apiKey = resolveOkpApiKey('', params.path);

if nargin < 3 || isempty(jobId)
    jobId = readJobIdFromMeta(fullfile(params.path,'data','OKP_job.txt'));
end

base    = 'https://predictor.openkinetics.org/api/v1';
opts    = weboptions('HeaderFields', {'Authorization', ['Bearer ' apiKey]}, ...
                     'Timeout', 60, 'ContentType', 'json');

while true
    statusResp = webread([base '/status/' jobId '/'], opts);
    state = lower(string(statusResp.status));

    switch state
        case "completed"
            csvText = webread([base '/result/' jobId '/'], ...
                weboptions('HeaderFields', {'Authorization', ['Bearer ' apiKey]}, ...
                           'Timeout', 120, 'ContentType', 'text'));
            fID = fopen(outFile,'w');
            fwrite(fID, csvText);
            fclose(fID);
            fprintf('OKP job %s completed; result stored at %s\n', jobId, outFile);
            kcatList = parseOkpOutput(model, outFile);
            done = true;
            return

        case "failed"
            error('OKP job %s failed. Check https://predictor.openkinetics.org/ for details.', jobId);

        otherwise % pending / running / queued
            pct = '';
            if isfield(statusResp,'progress') && isfield(statusResp.progress,'predictionsMade') ...
                    && isfield(statusResp.progress,'predictionsTotal') && statusResp.progress.predictionsTotal > 0
                pct = sprintf(' (%d/%d predictions)', statusResp.progress.predictionsMade, statusResp.progress.predictionsTotal);
            end
            if wait
                fprintf('OKP job %s status: %s%s. Waiting %d s...\n', jobId, char(statusResp.status), pct, pollInterval);
                pause(pollInterval);
            else
                fprintf(['OKP job %s not finished (status: %s%s).\n' ...
                         'Try again later with fetchOpenKineticsPredictor, ' ...
                         'or check https://predictor.openkinetics.org/.\n'], ...
                         jobId, char(statusResp.status), pct);
                done = false;
                kcatList = [];
                return
            end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function jobId = readJobIdFromMeta(metaFile)
if ~exist(metaFile,'file')
    error(['No jobId provided and no metadata file at %s. Run ' ...
           'submitOpenKineticsPredictor first, or pass jobId explicitly.'], metaFile);
end
lines = strsplit(fileread(metaFile), newline);
jobId = '';
for i=1:numel(lines)
    tok = regexp(strtrim(lines{i}), '^jobId:\s*(\S+)$', 'tokens', 'once');
    if ~isempty(tok); jobId = tok{1}; break; end
end
if isempty(jobId)
    error('Could not read a jobId from %s.', metaFile);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kcatList = parseOkpOutput(model, outFile)
% Parse the OKP result CSV into a kcatList. (Former
% readOpenKineticsPredictorOutput logic.) Result columns:
% kcat (1/s), Source kcat, Extra Info kcat, Protein Sequence, Substrate
% [, mean/max similarity ...]. Only columns 1,2,4,5 are used.

raw = readcell(outFile, 'Delimiter', ',', 'TextType', 'char');
if size(raw,2) < 5
    error('OKP output file does not have the expected >=5 columns: %s', outFile)
end
raw = raw(2:end, :);  % drop header

kcatsCol      = raw(:,1);
sourceKcatCol = raw(:,2);
sequences     = raw(:,4);
smiles        = raw(:,5);

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

valid = ~isnan(kcats);
kcats         = kcats(valid);
sequences     = sequences(valid);
smiles        = smiles(valid);
sourceKcatCol = sourceKcatCol(valid);

if isempty(kcats)
    error('OKP file does not contain any numeric kcat values.')
end
if ~isfield(model,'metSmiles')
    error('model.metSmiles is required to map SMILES back to metabolites.')
end

if model.ec.geckoLight
    origRxns = extractAfter(model.ec.rxns, 4);
else
    origRxns = model.ec.rxns;
end
[foundEcRxn, ecRxnToModelRxn] = ismember(origRxns, model.rxns);
if ~all(foundEcRxn)
    error('Not all entries in model.ec.rxns could be matched to model.rxns.')
end

rxnsOut = cell(0,1); genesOut = cell(0,1); substratesOut = cell(0,1);
kcatsOut = zeros(0,1); sourceOut = cell(0,1); unmatched = 0;

for i = 1:numel(sequences)
    proteinIdxs = find(strcmp(model.ec.sequence, sequences{i}));
    metIdxs     = find(strcmp(model.metSmiles, smiles{i}));
    if isempty(proteinIdxs) || isempty(metIdxs)
        unmatched = unmatched + 1; continue
    end
    modelRxnsUsingSub = any(model.S(metIdxs, :) < 0, 1);
    ecRxnUsesSub      = modelRxnsUsingSub(ecRxnToModelRxn).';
    ecRxnCatalyzed    = any(model.ec.rxnEnzMat(:, proteinIdxs) ~= 0, 2);
    candidateEcRxns   = find(ecRxnUsesSub & ecRxnCatalyzed);
    if isempty(candidateEcRxns)
        unmatched = unmatched + 1; continue
    end
    for k = 1:numel(candidateEcRxns)
        ecIdx           = candidateEcRxns(k);
        catalysts       = find(model.ec.rxnEnzMat(ecIdx, :) ~= 0);
        matchingProtein = intersect(proteinIdxs, catalysts, 'stable');
        if ~isempty(matchingProtein)
            geneName = model.ec.genes{matchingProtein(1)};
        else
            geneName = '';
        end
        origRxnIdx       = ecRxnToModelRxn(ecIdx);
        metIsSub         = model.S(metIdxs, origRxnIdx) < 0;
        substrateMetIdxs = metIdxs(metIsSub);
        if ~isempty(substrateMetIdxs)
            subName = model.metNames{substrateMetIdxs(1)};
        else
            subName = '';
        end
        rxnsOut{end+1,1}       = model.ec.rxns{ecIdx};       %#ok<AGROW>
        genesOut{end+1,1}      = geneName;                   %#ok<AGROW>
        substratesOut{end+1,1} = subName;                    %#ok<AGROW>
        kcatsOut(end+1,1)      = kcats(i);                   %#ok<AGROW>
        sourceOut{end+1,1}     = sourceKcatCol{i};           %#ok<AGROW>
    end
end

if unmatched > 0
    warning('%d of %d OKP entries could not be matched to an ec.rxn.', unmatched, numel(sequences));
end
if isempty(rxnsOut)
    error('No OKP entries could be mapped to model.ec.rxns.')
end

% Per-entry provenance from the 'Source kcat' column. OKP returns
% predicted values ('Prediction from CataPro' -> 'CataPro') as well as
% experimental values from databases ('BRENDA', 'Sabio-RK', 'UniProt'),
% which are kept verbatim. This is what selectKcatValue writes into
% model.ec.source.
kcatSource = regexprep(sourceOut, '^Prediction from\s+', '');

% Scalar source label: the most frequent provenance among the entries
% (e.g. 'CataPro' for a typical run). selectKcatValue only uses this as a
% fallback when kcatSource is absent, which it is not here.
uniqueSources = unique(kcatSource);
srcCounts     = cellfun(@(s) sum(strcmp(kcatSource, s)), uniqueSources);
[~, mostIdx]  = max(srcCounts);

kcatList.source     = uniqueSources{mostIdx};
kcatList.rxns       = rxnsOut;
kcatList.genes      = genesOut;
kcatList.substrates = substratesOut;
kcatList.kcats      = kcatsOut;
kcatList.kcatSource = kcatSource;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = resolveOkpApiKey(argKey, basePath)
if ~isempty(argKey)
    key = strtrim(char(argKey)); return
end
envKey = getenv('OKP_API_KEY');
if ~isempty(envKey)
    key = strtrim(envKey); return
end
keyFile = fullfile(basePath,'data','okpApiKey.txt');
if exist(keyFile,'file')
    key = strtrim(fileread(keyFile)); return
end
error(['No OpenKineticsPredictor API key found. Set the OKP_API_KEY ' ...
       'environment variable or place it in %s. Generate a key at ' ...
       'https://predictor.openkinetics.org/.'], keyFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = normalizeCellStrings(c)
for i = 1:numel(c)
    v = c{i};
    if ismissing(v)
        c{i} = '';
    elseif isstring(v)
        c{i} = char(v);
    elseif isnumeric(v)
        if isnan(v); c{i} = ''; else; c{i} = num2str(v); end
    elseif ~ischar(v)
        c{i} = char(string(v));
    end
end
end
