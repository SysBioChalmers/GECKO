function jobId = submitOpenKineticsPredictor(model, varargin)
% submitOpenKineticsPredictor  Submit a kcat prediction job to OpenKineticsPredictor.
%
% Builds the OpenKineticsPredictor input (protein sequences + single-substrate
% SMILES, via writeOpenKineticsPredictorInput) and submits a kcat prediction
% job to the OKP REST API (https://predictor.openkinetics.org/). The returned
% job id is also stored in data/OKP_job.txt, so fetchOpenKineticsPredictor can
% pick it up later.
%
% Predictors available via OKP: CataPro, CatPred, DLKcat, EITLEM, KinForm-H,
% KinForm-L, UniKP (see GET /api/v1/methods/).
%
% Obtaining an API key:
%   The OKP API requires a personal Bearer key (looks like 'ak_...'). Keys
%   are free and need no registration:
%     1. Open https://predictor.openkinetics.org/api-docs in a browser.
%     2. In the "API key generator" section, click "Generate".
%     3. Copy the key shown (it is displayed only once). If you lose it,
%        revoke it on the same page and generate a new one.
%   Keys are tied to your IP and carry a daily prediction quota that resets
%   at midnight UTC.
%
%   The key is a secret; do NOT put it in the model adapter (which is
%   shared/committed). Provide it in one of three ways, checked in this
%   order: the apiKey argument below; the OKP_API_KEY environment variable
%   (setenv('OKP_API_KEY',...)); or a plain-text file data/okpApiKey.txt in
%   the model's data folder (this name is git-ignored by GECKO).
%
% Parameters
% ----------
% model : struct
%     ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% ecRxns : logical
%     logical vector indicating which reactions to include (default: all
%     reactions).
% overwrite : logical
%     rebuild data/OKP.csv even if it exists (default false; an existing
%     file is reused).
% apiKey : char
%     OKP API key (default: resolved from OKP_API_KEY or data/okpApiKey.txt,
%     see above).
% method : char
%     predictor to use (default: params.okp.method, or 'CataPro').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% jobId : char
%     the OKP job identifier (also written to data/OKP_job.txt).
%
% Examples
% --------
%     jobId = submitOpenKineticsPredictor(model);
%     jobId = submitOpenKineticsPredictor(model, 'overwrite', true, 'apiKey', 'ak_...', 'method', 'DLKcat');
%
% See also
% --------
% fetchOpenKineticsPredictor, writeOpenKineticsPredictorInput

p = parseGECKOargs(varargin, { ...
    'ecRxns',       []; ...
    'overwrite',    []; ...
    'apiKey',       []; ...
    'method',       []; ...
    'modelAdapter', []});
ecRxns       = p.ecRxns;
overwrite    = p.overwrite;
apiKey       = p.apiKey;
method       = p.method;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either provide a modelAdapter or set the default in ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(overwrite)
    overwrite = false;
end
if isempty(ecRxns)
    ecRxns = true(numel(model.ec.rxns),1);
elseif ~islogical(ecRxns) || numel(ecRxns) ~= numel(model.ec.rxns)
    error('ecRxns should be a logical vector with length equal to model.ec.rxns')
end

apiKey = resolveOkpApiKey(apiKey, params.path);

% Resolve predictor method and request options: defaults -> params.okp.* ->
% the method argument (highest priority, applied last).
okp = struct('method','CataPro', ...
             'handleLongSequences','truncate', ...
             'includeSimilarityColumns',true,'canonicalizeSubstrates',true);
if isfield(params,'okp')
    fn = fieldnames(params.okp);
    for i=1:numel(fn); okp.(fn{i}) = params.okp.(fn{i}); end
end
if ~isempty(method)
    okp.method = method;
end

%% Build (or reuse) the input CSV --- the same function, writing the same
% file, as the file-based workflow (writeOpenKineticsPredictorInput /
% readOpenKineticsPredictorOutput), rather than a separate reimplementation.
filename = fullfile(params.path,'data','OKP.csv');
if exist(filename,'file') && ~overwrite
    fprintf('Using existing %s (set overwrite=true to rebuild).\n', filename);
else
    writeOpenKineticsPredictorInput(model, 'ecRxns', ecRxns, 'modelAdapter', modelAdapter, ...
        'onlyWithSmiles', true, 'overwrite', true);
end

%% Submit via the OKP REST API
% GECKO only ever predicts kcat.
url = 'https://predictor.openkinetics.org/api/v1/submit/';

targetsJson = '["kcat"]';
methodsJson = ['{"kcat":"' okp.method '"}'];

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
fprintf(fID,'submittedAt: %s\n', char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
fclose(fID);

fprintf('Submitted OKP job %s (method: %s). Metadata stored at %s\n', jobId, okp.method, metaFile);
fprintf('Check progress later with fetchOpenKineticsPredictor, or at https://predictor.openkinetics.org/\n');
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
