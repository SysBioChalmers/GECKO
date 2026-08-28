function [done, kcatList] = fetchOpenKineticsPredictor(model, varargin)
% fetchOpenKineticsPredictor  Fetch and parse an OpenKineticsPredictor job.
%
% Checks the status of an OpenKineticsPredictor job and, when it has
% finished, downloads the result to data/OKP_output.csv and parses it (via
% readOpenKineticsPredictorOutput) into a kcatList consumable by
% selectKcatValue.
%
% The API key (only needed when useStored=false) is resolved from the
% OKP_API_KEY environment variable or data/okpApiKey.txt; see
% submitOpenKineticsPredictor for how to obtain and store a key.
%
% Parameters
% ----------
% model : struct
%     ecModel in GECKO 3 format (with ecModel.ec structure; model.metSmiles
%     is required to map results back).
%
% Name-Value Arguments
% --------------------
% useStored : logical
%     if true, skip the API entirely and parse the already-downloaded
%     data/OKP_output.csv. If false, query the API (default false).
% jobId : char
%     OKP job id (default: read from data/OKP_job.txt, written by
%     submitOpenKineticsPredictor).
% wait : logical
%     if true, poll until the job finishes; if false, report once and
%     return (default false).
% pollInterval : double
%     seconds between polls when wait=true (default 30).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% done : logical
%     true if a result was obtained (downloaded or loaded from disk).
% kcatList : struct
%     structure with kcat values (empty when not done); see
%     readOpenKineticsPredictorOutput for its fields (source, rxns, genes,
%     substrates, kcats, kcatSource).
%
% Examples
% --------
%     [done, kcatList] = fetchOpenKineticsPredictor(model);                        % check once
%     [done, kcatList] = fetchOpenKineticsPredictor(model, 'wait', true);          % poll
%     [done, kcatList] = fetchOpenKineticsPredictor(model, 'useStored', true);     % parse stored file
%
% See also
% --------
% submitOpenKineticsPredictor, readOpenKineticsPredictorOutput

p = parseGECKOargs(varargin, { ...
    'useStored',    []; ...
    'jobId',        []; ...
    'wait',         []; ...
    'pollInterval', []; ...
    'modelAdapter', []});
useStored    = p.useStored;
jobId        = p.jobId;
wait         = p.wait;
pollInterval = p.pollInterval;
modelAdapter = p.modelAdapter;

if isempty(useStored); useStored = false; end
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
if isempty(wait); wait = false; end
if isempty(pollInterval); pollInterval = 30; end

outFile = fullfile(params.path,'data','OKP_output.csv');

%% Stored-result mode: parse the saved file, no API call
if useStored
    if ~exist(outFile,'file')
        error('useStored=true but no stored result found at %s. Run fetchOpenKineticsPredictor with useStored=false first.', outFile);
    end
    kcatList = readOpenKineticsPredictorOutput(model, 'outFile', outFile, 'modelAdapter', modelAdapter);
    done = true;
    return
end

%% API mode
apiKey = resolveOkpApiKey('', params.path);

if isempty(jobId)
    jobId = readJobIdFromMeta(fullfile(params.path,'data','OKP_job.txt'));
end

base = 'https://predictor.openkinetics.org/api/v1';
opts = weboptions('HeaderFields', {'Authorization', ['Bearer ' apiKey]}, ...
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
            kcatList = readOpenKineticsPredictorOutput(model, 'outFile', outFile, 'modelAdapter', modelAdapter);
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
