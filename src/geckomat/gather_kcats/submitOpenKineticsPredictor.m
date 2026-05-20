function jobId = submitOpenKineticsPredictor(model, ecRxns, modelAdapter, method, apiKey, overwrite)
% submitOpenKineticsPredictor
%   Builds the OpenKineticsPredictor input (protein sequences + single-
%   substrate SMILES) and submits a prediction job directly to the OKP
%   REST API (https://predictor.openkinetics.org/), replacing the manual
%   upload step. The returned job id is also stored in data/OKP_job.txt,
%   so fetchOpenKineticsPredictor can pick it up later.
%
%   Predictors available via OKP: CataPro, CatPred, DLKcat, EITLEM,
%   KinForm-H, KinForm-L, UniKP (see GET /api/v1/methods/).
%
% Input:
%   model           ecModel in GECKO 3 format (with ecModel.ec structure)
%   ecRxns          logical vector indicating which reactions to include
%                   (Optional, default: all reactions)
%   modelAdapter    loaded model adapter (Optional, uses default if empty)
%   method          predictor to use (Optional). Resolved as: this arg ->
%                   params.okp.method -> 'CataPro'.
%   apiKey          OKP API key (Optional). Resolved as: this arg ->
%                   environment variable OKP_API_KEY -> data/okpApiKey.txt.
%                   Generate a key at https://predictor.openkinetics.org/.
%   overwrite       logical, rebuild data/OKP.csv even if it exists
%                   (Optional, default: false; an existing file is reused).
%
% Output:
%   jobId           the OKP job identifier (also written to data/OKP_job.txt)
%
% Usage:
%   jobId = submitOpenKineticsPredictor(model);
%   jobId = submitOpenKineticsPredictor(model, ecRxns, modelAdapter, 'DLKcat');

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

% Resolve OKP settings: function arg -> params.okp.* -> built-in default
okp = struct('method','CataPro','targets',{{'kcat'}}, ...
             'handleLongSequences','truncate', ...
             'includeSimilarityColumns',true,'canonicalizeSubstrates',true);
if isfield(params,'okp')
    fn = fieldnames(params.okp);
    for i=1:numel(fn); okp.(fn{i}) = params.okp.(fn{i}); end
end
if nargin >= 4 && ~isempty(method)
    okp.method = method;
end

if nargin < 5; apiKey = ''; end
apiKey = resolveOkpApiKey(apiKey, params.path);

if nargin<6 || isempty(overwrite)
    overwrite = false;
end

%% Build (or reuse) the input CSV
filename = fullfile(params.path,'data','OKP.csv');
if exist(filename,'file') && ~overwrite
    fprintf('Using existing %s (set overwrite=true to rebuild).\n', filename);
else
    buildOkpInput(model, ecRxns, params, geckoPath, filename);
end

%% Submit via the OKP REST API
url = 'https://predictor.openkinetics.org/api/v1/submit/';

targetsJson = ['["' strjoin(cellstr(okp.targets),'","') '"]'];
methodsJson = ['{"' strjoin(cellstr(okp.targets),['":"' okp.method '","']) '":"' okp.method '"}'];

import matlab.net.http.*
import matlab.net.http.io.*
provider = MultipartFormProvider( ...
    'file',                     FileProvider(filename), ...
    'targets',                  targetsJson, ...
    'methods',                  methodsJson, ...
    'handleLongSequences',      okp.handleLongSequences, ...
    'includeSimilarityColumns', boolToStr(okp.includeSimilarityColumns), ...
    'canonicalizeSubstrates',   boolToStr(okp.canonicalizeSubstrates));
header = matlab.net.http.HeaderField('Authorization', ['Bearer ' apiKey]);
req    = RequestMessage('POST', header, provider);
resp   = req.send(url);

data = resp.Body.Data;
if resp.StatusCode ~= matlab.net.http.StatusCode.OK && ...
        resp.StatusCode ~= matlab.net.http.StatusCode.Created
    error('OKP submit failed (HTTP %d): %s', double(resp.StatusCode), okpErrorText(data));
end
if ~isstruct(data) || ~isfield(data,'jobId')
    error('OKP submit returned an unexpected response (no jobId).');
end
jobId = data.jobId;

%% Persist job metadata (plain text, one field per line)
metaFile = fullfile(params.path,'data','OKP_job.txt');
fID = fopen(metaFile,'w');
fprintf(fID,'jobId: %s\n', jobId);
fprintf(fID,'method: %s\n', okp.method);
fprintf(fID,'targets: %s\n', strjoin(cellstr(okp.targets),','));
fprintf(fID,'submittedAt: %s\n', char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
fclose(fID);

fprintf('Submitted OKP job %s (method: %s). Metadata stored at %s\n', jobId, okp.method, metaFile);
fprintf('Check progress later with fetchOpenKineticsPredictor, or at https://predictor.openkinetics.org/\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buildOkpInput(model, ecRxns, params, geckoPath, filename)
% Build data/OKP.csv with header 'Protein Sequence,Substrate', one row
% per (enzyme sequence, substrate SMILES) pair for the selected reactions.
% (Former writeOpenKineticsPredictorInput logic.)

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

% Metabolites to ignore (reuses DLKcat files)
if exist(fullfile(params.path,'data','DLKcatIgnoreMets.tsv'),'file')
    ignoreFile = fullfile(params.path,'data','DLKcatIgnoreMets.tsv');
else
    ignoreFile = fullfile(geckoPath,'databases','DLKcatIgnoreMets.tsv');
end
fID = fopen(ignoreFile);
fileData = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[ignoreMets, ignoreSmiles] = deal(fileData{[1,2]});

metsNorm        = lower(regexprep(model.metNames,'[^0-9a-zA-Z]+',''));
ignoreMetsNorm  = lower(regexprep(ignoreMets,'[^0-9a-zA-Z]+',''));
ignoreSmiles(cellfun(@isempty,ignoreSmiles)) = [];

ignoreMetsIdx = ismember(metsNorm, ignoreMetsNorm);
if isfield(model,'metSmiles')
    ignoreMetsIdx = ignoreMetsIdx | ismember(model.metSmiles, ignoreSmiles);
end
ignoreMetsIdx = ignoreMetsIdx | startsWith(model.mets,'prot_');

reducedS = model.S;
reducedS(ignoreMetsIdx,:) = 0;

% Remove currency metabolite pairs (reuses DLKcat files)
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
    rxnsWithRemainingSubstrates = any(tempRedS(:,pairRxns) < 0,1);
    reducedS([find(subs);find(prod)],intersect(pairRxns,pairRxns(rxnsWithRemainingSubstrates))) = 0;
end

% Extract substrates for selected reactions only
clearedS = reducedS(:, origRxnIdxs);
[substrateIdxs, reactionIdxs] = find(clearedS < 0);
[proteinIdxs, rxnEnzIdxs] = find(model.ec.rxnEnzMat(reactionIdxs,:)');

sequences = model.ec.sequence(proteinIdxs);
if isfield(model,'metSmiles')
    smiles = model.metSmiles(substrateIdxs(rxnEnzIdxs));
else
    smiles = cell(size(sequences));
end

% Only include entries with both sequence and SMILES
valid = ~cellfun(@isempty, sequences) & ~cellfun(@isempty, smiles);
sequences = sequences(valid);
smiles    = smiles(valid);
if isempty(sequences)
    error('No (sequence, SMILES) pairs to submit. Ensure model.metSmiles is populated.')
end

outTable = [sequences(:)'; smiles(:)'];
fID = fopen(filename,'w');
fprintf(fID,'Protein Sequence,Substrate\n');
fprintf(fID,'%s,%s\n', outTable{:});
fclose(fID);
fprintf('OpenKineticsPredictor input stored at %s\n', filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = resolveOkpApiKey(argKey, basePath)
% Resolve the OKP API key: function arg -> OKP_API_KEY env -> data/okpApiKey.txt
if ~isempty(argKey)
    key = strtrim(char(argKey));
    return
end
envKey = getenv('OKP_API_KEY');
if ~isempty(envKey)
    key = strtrim(envKey);
    return
end
keyFile = fullfile(basePath,'data','okpApiKey.txt');
if exist(keyFile,'file')
    key = strtrim(fileread(keyFile));
    return
end
error(['No OpenKineticsPredictor API key found. Provide it as the apiKey ' ...
       'argument, set the OKP_API_KEY environment variable, or place it in ' ...
       '%s. Generate a key at https://predictor.openkinetics.org/.'], keyFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = boolToStr(b)
if islogical(b) || isnumeric(b)
    if b; s = 'true'; else; s = 'false'; end
else
    s = char(b);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = okpErrorText(data)
if isstruct(data) && isfield(data,'error')
    txt = data.error;
elseif ischar(data) || isstring(data)
    txt = char(data);
else
    txt = 'unknown error';
end
end
