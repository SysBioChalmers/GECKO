function [kcatList, done] = kcatsFromOKP(model, varargin)
% kcatsFromOKP  Gather kcats from OpenKineticsPredictor and keep them as a source file.
%
% Predicts kcats with one OpenKineticsPredictor (OKP) method and writes
% them to data/kcats/okp_<method>.tsv in the project folder (see
% saveKcatSource). OKP is a web service with a job queue, so this function
% can be called repeatedly:
%
% - without a job for this method, it submits one
%   (submitOpenKineticsPredictor);
% - with a job, it fetches the result (fetchOpenKineticsPredictor). When
%   the job is done the file is written; otherwise it reports the job's
%   status and writes nothing.
%
% Each method keeps its own job in data/OKP_job_<method>.txt and its own
% result in data/OKP_output_<method>.csv, so jobs of several methods can
% run at the same time. The job file records which job the stored result
% came from; useStored only uses a result of the current job. An API key
% is needed (see resolveOkpApiKey).
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure), with
%     metabolite SMILES (findMetSmiles) and protein sequences.
%
% Name-Value Arguments
% --------------------
% method : char
%     OKP method, e.g. 'CataPro', 'CatPred', 'CatRange', 'DLKcat',
%     'EITLEM', 'KinForm-H', 'KinForm-L', 'OmniESI', 'TurNuP', 'UniKP'
%     (default: params.okp.method of the model adapter, or 'CataPro').
% wait : logical
%     wait until the job is done (default false).
% resubmit : logical
%     submit a new job even if one exists for this method; its earlier
%     result is removed (default false).
% useStored : logical
%     read the stored result of this method's job, without contacting OKP
%     (default false).
% ecRxns : logical
%     for which model.ec.rxns kcats are predicted (default: all). Only
%     used when a job is submitted.
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'okp_<method>', in
%     lower case, with other characters than letters and digits as '_').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values, as returned by readOpenKineticsPredictorOutput;
%     empty while the job is not done.
% done : logical
%     whether the result was available.
%
% Examples
% --------
%     kcatsFromOKP(ecModel, 'method', 'CataPro');          % submits
%     kcatsFromOKP(ecModel, 'method', 'CataPro');          % later: fetches and writes
%     kcatsFromOKP(ecModel, 'method', 'UniKP', 'wait', true);
%
% See also
% --------
% submitOpenKineticsPredictor, fetchOpenKineticsPredictor, kcatsFromBrenda, saveKcatSource

checkNameValue(varargin, {'method','wait','resubmit','useStored','ecRxns','save','saveAs','modelAdapter'}, 'kcatsFromOKP');
p = parseGECKOargs(varargin, { ...
    'method',       []; ...
    'wait',         false; ...
    'resubmit',     false; ...
    'useStored',    false; ...
    'ecRxns',       []; ...
    'save',         true; ...
    'saveAs',       []; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
method = char(p.method);
if isempty(method)
    if isfield(params, 'okp') && isfield(params.okp, 'method') && ~isempty(params.okp.method)
        method = char(params.okp.method);
    else
        method = 'CataPro';
    end
end
tag = lower(regexprep(method, '[^A-Za-z0-9]+', '_'));
dataDir = fullfile(params.path, 'data');
jobFile = fullfile(dataDir, ['OKP_job_' tag '.txt']);
outFile = fullfile(dataDir, ['OKP_output_' tag '.csv']);
name = char(p.saveAs);
if isempty(name)
    name = ['okp_' tag];
end
job = readJobFile(jobFile);
submitting = ~p.useStored && (p.resubmit || ~isfield(job, 'jobId'));
if ~isempty(p.ecRxns) && ~submitting
    error('kcatsFromOKP:ecRxns', ['ecRxns selects the reactions of a new job, but this call ' ...
        'submits none (a job for %s exists; use ''resubmit'', true for a new one).'], method)
end
if p.save
    file = kcatSourceFile(name, [], modelAdapter); % fails before submitting, not after
    if submitting
        warnSubset(p.ecRxns, numel(model.ec.rxns), file);
    end
end

if p.useStored
    if ~isfield(job, 'jobId')
        adoptSharedResult(model, modelAdapter, dataDir, method, jobFile, outFile);
        job = readJobFile(jobFile);
    end
    if ~isfile(outFile) || ~isfield(job, 'resultOf') || ~strcmp(job.resultOf, job.jobId)
        error('kcatsFromOKP:noStored', ['No stored result of the current %s job (%s). Call ' ...
            'kcatsFromOKP without ''useStored'' to fetch it.'], method, jobFile)
    end
    kcatList = readOpenKineticsPredictorOutput(model, 'outFile', outFile, 'modelAdapter', modelAdapter);
    done = true;
else
    if submitting
        if isfile(outFile)
            delete(outFile);
        end
        jobId = submitOpenKineticsPredictor(model, 'method', method, 'ecRxns', p.ecRxns, ...
            'overwrite', true, 'modelAdapter', modelAdapter);
        job = struct('jobId', jobId, 'method', method, 'submittedAt', ...
            char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
        if ~isempty(p.ecRxns)
            job.ecRxns = sprintf('%d', countRxns(p.ecRxns));
        end
        writeJobFile(jobFile, job);
    end
    [done, kcatList] = fetchOpenKineticsPredictor(model, 'jobId', job.jobId, 'wait', p.wait, ...
        'outFile', outFile, 'modelAdapter', modelAdapter);
    if done
        job.resultOf = job.jobId;
        writeJobFile(jobFile, job);
    end
end

if ~done
    fprintf('Call kcatsFromOKP(model, ''method'', ''%s'') again later to collect the result.\n', method);
    return
end
if p.save
    parameters = struct('method', method);
    if isfield(job, 'ecRxns')
        parameters.ec_rxns = str2double(job.ecRxns);
    end
    saveKcatSource(model, kcatList, name, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromOKP', 'parameters', parameters, 'data', ['OKP job ' job.jobId]);
    fprintf('OKP (%s) kcats for %d reactions written to %s\n', method, ...
        numel(unique(kcatList.rxns(kcatList.kcats > 0))), file);
end
end

function adoptSharedResult(model, modelAdapter, dataDir, method, jobFile, outFile)
% A result that fetchOpenKineticsPredictor stored by itself, in
% data/OKP_output.csv with data/OKP_job.txt, is taken over for this method
% only when its job is for this method and its predictions are labelled
% with this method.
shared = fullfile(dataDir, 'OKP_output.csv');
sharedJob = readJobFile(fullfile(dataDir, 'OKP_job.txt'));
ok = isfile(shared) && isfield(sharedJob, 'jobId') && isfield(sharedJob, 'method') && ...
    strcmpi(sharedJob.method, method);
if ok
    list = readOpenKineticsPredictorOutput(model, 'outFile', shared, 'modelAdapter', modelAdapter);
    predicted = list.kcatSource(startsWith(list.kcatSource, 'OKP-') & ...
        ~ismember(lower(list.kcatSource), {'okp-brenda','okp-sabio-rk','okp-uniprot'}));
    ok = ~isempty(predicted) && all(strcmpi(predicted, ['OKP-' method]));
end
if ~ok
    error('kcatsFromOKP:noStored', ['No stored OKP result for %s (%s). Call kcatsFromOKP ' ...
        'without ''useStored'' to submit or fetch a job.'], method, outFile)
end
copyfile(shared, outFile);
sharedJob.method = method;
sharedJob.resultOf = sharedJob.jobId;
writeJobFile(jobFile, sharedJob);
end

function job = readJobFile(file)
% Fields of an OKP job file ('key: value' per line), as a struct.
job = struct();
if ~isfile(file)
    return
end
lines = strsplit(fileread(file), newline);
for i = 1:numel(lines)
    tok = regexp(lines{i}, '^([A-Za-z]\w*):\s*(.*?)\s*$', 'tokens', 'once');
    if ~isempty(tok)
        job.(tok{1}) = tok{2};
    end
end
end

function writeJobFile(file, job)
fid = fopen(file, 'w');
for f = {'jobId','method','submittedAt','ecRxns','resultOf'}
    if isfield(job, f{1})
        fprintf(fid, '%s: %s\n', f{1}, job.(f{1}));
    end
end
fclose(fid);
end
