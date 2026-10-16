function report = updateKcatSources(model, varargin)
% updateKcatSources  Gather the kcat sources that the configuration needs.
%
% Lists the kcat source files that [kcat_assignment] in model_adapter.toml
% uses (see collectKcats), and gathers those that are missing or outdated
% with kcatsFromBrenda, kcatsFromDLKcat and kcatsFromOKP:
%
% - a file is outdated when its header records other settings than
%   [kcat_sources.<name>] asks for (the preset sets these for brenda), or,
%   for brenda, another BRENDA release than GECKO's BRENDA files;
% - brenda is matched again; dlkcat reads data/DLKcat.tsv if it exists and
%   runs DLKcat (Docker) otherwise; okp_<method> submits or fetches a job,
%   and reports it as pending while the job runs;
% - custom.tsv and other files it cannot gather are left as they are.
%
% Requirements are checked first: Docker for a DLKcat run and an
% OpenKineticsPredictor API key for OKP. Calling it again is safe: files
% that are current are not gathered again.
%
% Per-source settings are read from model_adapter.toml, for instance:
%
%     [kcat_sources.brenda]
%     aggregation = "median"
%
%     [kcat_sources.okp_catapro]
%     method = "CataPro"
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% dryRun : logical
%     only report what would be gathered (default false).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% report : table
%     one row per source: source, status (current, gathered, pending,
%     outdated, missing, failed or skipped) and note.
%
% Examples
% --------
%     updateKcatSources(ecModel, 'dryRun', true)
%     updateKcatSources(ecModel);
%
% See also
% --------
% collectKcats, assignKcats, kcatsFromBrenda, kcatsFromDLKcat, kcatsFromOKP

checkNameValue(varargin, {'dryRun','modelAdapter'}, 'updateKcatSources');
p = parseGECKOargs(varargin, {'dryRun', false; 'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
config = kcatAssignmentConfig(params);
kcatDir = fullfile(params.path, 'data', 'kcats');

names = {};
for s = 1:numel(config.stages)
    for t = 1:numel(config.stages(s).sources)
        token = config.stages(s).sources{t};
        if startsWith(token, 'database_')
            token = 'brenda';
        end
        if ~any(strcmp(names, token))
            names{end+1} = token; %#ok<AGROW>
        end
    end
end

status = cell(numel(names), 1);
note = cell(numel(names), 1);
for k = 1:numel(names)
    name = names{k};
    file = fullfile(kcatDir, [name '.tsv']);
    wanted = struct();
    if isfield(config.sources, name)
        wanted = config.sources.(name);
    end
    [status{k}, note{k}] = check(name, file, wanted);
    if any(strcmp(status{k}, {'missing','outdated'})) && ~p.dryRun
        [status{k}, note{k}] = gather(model, name, wanted, modelAdapter);
    end
end
report = table(names(:), status, note, 'VariableNames', {'source','status','note'});
disp(report)

    function [st, nt] = check(name, file, wanted)
        st = 'current'; nt = '';
        if strcmp(name, 'custom')
            st = 'skipped'; nt = 'curated by hand';
            return
        elseif ~any(strcmp(name, {'brenda','dlkcat'})) && ~startsWith(name, 'okp_')
            st = 'skipped'; nt = 'not a source GECKO gathers';
            return
        end
        if ~isfile(file)
            st = 'missing';
            return
        end
        [~, meta] = readKcatSource(file);
        recorded = struct();
        if isfield(meta, 'parameters') && ~isempty(meta.parameters)
            recorded = jsondecode(meta.parameters);
        end
        f = fieldnames(wanted);
        for i = 1:numel(f)
            if ~isfield(recorded, f{i}) || ~isequal(normalizeValue(recorded.(f{i})), normalizeValue(wanted.(f{i})))
                st = 'outdated';
                nt = sprintf('%s differs from model_adapter.toml', f{i});
                return
            end
        end
        if strcmp(name, 'brenda')
            release = brendaRelease(modelAdapter);
            if isfield(meta, 'data') && ~isempty(release) && ~strcmp(meta.data, release)
                st = 'outdated';
                nt = sprintf('made from %s; GECKO has %s', meta.data, release);
            end
        end
    end
end

function [st, nt] = gather(model, name, wanted, modelAdapter)
st = 'gathered'; nt = '';
try
    if strcmp(name, 'brenda')
        args = {};
        if isfield(wanted, 'aggregation'), args = [args, {'aggregation', wanted.aggregation}]; end
        if isfield(wanted, 'force_wildcard_level'), args = [args, {'forceWClvl', wanted.force_wildcard_level}]; end
        kcatsFromBrenda(model, args{:}, 'modelAdapter', modelAdapter);
    elseif strcmp(name, 'dlkcat')
        run = ~isfile(fullfile(modelAdapter.params.path, 'data', 'DLKcat.tsv'));
        if run
            [status, ~] = system('docker --version');
            if status ~= 0
                st = 'failed'; nt = 'DLKcat needs Docker to run';
                return
            end
        end
        kcatsFromDLKcat(model, 'run', run, 'modelAdapter', modelAdapter);
    else % okp_<method>
        method = okpMethod(name, wanted);
        resolveOkpApiKey('', modelAdapter.params.path); % errors without a key
        [~, done] = kcatsFromOKP(model, 'method', method, 'saveAs', name, 'modelAdapter', modelAdapter);
        if ~done
            st = 'pending'; nt = sprintf('OKP job for %s still running; call again later', method);
        end
    end
catch err
    st = 'failed'; nt = err.message;
end
end

function method = okpMethod(name, wanted)
% The OKP method of okp_<method>: from [kcat_sources.okp_<method>] method,
% or the known method whose name gives this file name.
if isfield(wanted, 'method')
    method = char(wanted.method);
    return
end
known = {'CataPro','CatPred','CatRange','DLKcat','EITLEM','KinForm-H','KinForm-L', ...
    'OmniESI','OmniESI-O2DENet','TurNuP','UniKP'};
tags = cellfun(@(m) ['okp_' lower(regexprep(m, '[^A-Za-z0-9]+', '_'))], known, 'UniformOutput', false);
hit = strcmp(tags, name);
if ~any(hit)
    error('updateKcatSources:okpMethod', ['Set the OKP method of %s in model_adapter.toml: ' ...
        '[kcat_sources.%s] method = "..."'], name, name)
end
method = known{hit};
end

function v = normalizeValue(v)
if isstring(v)
    v = char(v);
end
if ischar(v)
    v = lower(v);
elseif islogical(v)
    v = double(v);
end
end

function release = brendaRelease(modelAdapter)
release = '';
fid = fopen(fullfile(modelAdapter.getBrendaDBFolder(), 'kcat.tsv'), 'r');
if fid < 0
    return
end
line = fgetl(fid);
fclose(fid);
tok = regexp(line, '^#\s*(BRENDA release [^ ]+)', 'tokens', 'once');
if ~isempty(tok)
    release = tok{1};
end
end
