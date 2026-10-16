function [decisions, report, ecModel] = collectKcats(model, varargin)
% collectKcats  Decide one kcat per reaction from the kcat source files.
%
% Combines the kcat source files in data/kcats (see kcatsFromBrenda,
% kcatsFromDLKcat, kcatsFromOKP and saveKcatSource) as [kcat_assignment] in
% model_adapter.toml sets, without changing the model given. assignKcats
% applies the decisions.
%
% The order is a preset or a list of stages:
%
%     [kcat_assignment]
%     preset = "gecko_median"      # or gecko_classic (the default)
%     isozymes = true              # fill from isozymes after the stages
%     standard = true              # standard kcat for the rest
%
%     # instead of a preset:
%     [[kcat_assignment.stage]]
%     action = "select"
%     sources = ["custom", "database_top", "dlkcat", "database_bottom"]
%     criteria = "max"             # among several values of one source
%     top_origin_limit = 6
%     bottom_origin_limit = 6
%     wildcard_limit = 3
%
% A select stage decides each reaction that is still open by the first
% listed source with a value (see mergeKcats). Sources are the names of
% source files (custom, dlkcat, okp_catapro, ...) and the BRENDA quality
% tiers database_exact, database_top and database_bottom, which read
% brenda.tsv. Both presets order custom, database_top, dlkcat and
% database_bottom; gecko_classic uses the BRENDA maximum and gecko_median
% the BRENDA median ([kcat_sources.brenda] aggregation, see
% updateKcatSources).
%
% After the stages:
% - custom kcats of data/customKcats.tsv are applied (applyCustomKcats),
%   when a stage lists custom but there is no data/kcats/custom.tsv;
% - with isozymes, reactions without a kcat take one from their isozymes
%   (getKcatAcrossIsozymes);
% - with standard, the remaining reactions get the standard kcat
%   (getStandardKcat).
%
% A source file a stage names that does not exist is reported and left out.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% writeReport : logical
%     write the decisions to output/<model id>_kcatDecisions.tsv (default
%     true).
%
% Returns
% -------
% decisions : table
%     one row per model.ec reaction: ecRxn, stage (the stage number, or
%     custom, isozymes, standard, or empty when no kcat was found), source
%     (as in model.ec.source), kcat, candidates (the number of values the
%     deciding source had for the reaction) and spreadLog10 (their range,
%     in log10 units).
% report : struct
%     - preset: the preset used, or '' with stages;
%     - missing: source files the stages name that do not exist;
%     - sources: per source file loaded, its unmatched, partial and invalid
%       rows (see loadKcatSource).
% ecModel : struct
%     the model with the decided kcats, before applyKcatConstraints; this
%     is what assignKcats returns, after applying the constraints.
%
% Examples
% --------
%     decisions = collectKcats(ecModel);
%     decisions(strcmp(decisions.stage, ''), :)   % reactions without a kcat
%
% See also
% --------
% assignKcats, updateKcatSources, loadKcatSource, mergeKcats, assignKcatValues

checkNameValue(varargin, {'modelAdapter','writeReport'}, 'collectKcats');
p = parseGECKOargs(varargin, {'modelAdapter', []; 'writeReport', true});
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

ecModel = model;
nEc = numel(ecModel.ec.rxns);
ecModel.ec.kcat(:) = 0;
ecModel.ec.source(:) = {''};
stage = repmat({''}, nEc, 1);
candidates = zeros(nEc, 1);
spread = nan(nEc, 1);

report.preset = config.preset;
report.missing = {};
report.sources = struct();
loaded = struct();
legacyCustom = false;

    function list = loadSource(name)
        % The kcat list of a source file, or [] when it does not exist.
        list = [];
        file = fullfile(kcatDir, [name '.tsv']);
        if isfile(file)
            [list, rep] = loadKcatSource(ecModel, file);
            report.sources.(name) = rep;
        elseif strcmp(name, 'custom') && isfile(fullfile(params.path, 'data', 'customKcats.tsv'))
            legacyCustom = true;
        else
            report.missing{end+1} = name;
            printOrange(sprintf(['No kcat source file %s; it is left out. ' ...
                'updateKcatSources gathers it.\n'], file));
        end
    end

for s = 1:numel(config.stages)
    st = config.stages(s);
    lists = {};
    priority = {};
    for t = 1:numel(st.sources)
        token = st.sources{t};
        name = sourceFile(token);
        if ~isfield(loaded, name)
            loaded.(name) = loadSource(name);
        end
        list = loaded.(name);
        if isempty(list)
            continue
        end
        if ~any(cellfun(@(l) isequal(l, list), lists))
            lists{end+1} = list; %#ok<AGROW>
        end
        if startsWith(token, 'database_')
            priority{end+1} = token; %#ok<AGROW>
        else
            % The labels of the file's rows, at the file's place in the order.
            labels = unique(cellfun(@normalizeSource, list.kcatSource, 'UniformOutput', false), 'stable');
            priority = [priority, labels(~ismember(labels, {'brenda','sabio_rk'}))']; %#ok<AGROW>
        end
    end
    if isempty(lists) || isempty(priority)
        continue
    end
    merged = mergeKcats(lists, priority, 'topOriginLimit', st.topOriginLimit, ...
        'bottomOriginLimit', st.bottomOriginLimit, 'wildcardLimit', st.wildcardLimit);
    if isempty(merged.rxns)
        continue
    end
    open = ecModel.ec.kcat == 0;
    [ecModel, idx] = assignKcatValues(ecModel, merged, 'criteria', st.criteria, 'overwrite', 'false');
    idx = idx(open(idx));
    stage(idx) = {sprintf('%d', s)};
    [~, where] = ismember(merged.rxns, ecModel.ec.rxns);
    valid = merged.kcats > 0;
    for i = idx(:)'
        values = merged.kcats(where == i & valid);
        candidates(i) = numel(values);
        spread(i) = log10(max(values)) - log10(min(values));
    end
end

if legacyCustom
    before = ecModel.ec.kcat;
    ecModel = applyCustomKcats(ecModel, 'modelAdapter', modelAdapter);
    changed = ecModel.ec.kcat ~= before;
    stage(changed) = {'custom'};
    candidates(changed) = 1;
    spread(changed) = 0;
end
if config.isozymes && ~ecModel.ec.geckoLight
    before = ecModel.ec.kcat == 0;
    ecModel = getKcatAcrossIsozymes(ecModel);
    stage(before & ecModel.ec.kcat > 0) = {'isozymes'};
end
if config.standard
    ecModel = getStandardKcat(ecModel, 'modelAdapter', modelAdapter);
    n = numel(ecModel.ec.rxns);
    stage(end+1:n) = {''};
    candidates(end+1:n) = 0;
    spread(end+1:n) = NaN;
    stage(strcmp(ecModel.ec.source, 'standard') & ecModel.ec.kcat > 0) = {'standard'};
end

decisions = table(ecModel.ec.rxns(:), stage(:), ecModel.ec.source(:), ecModel.ec.kcat(:), ...
    candidates(:), spread(:), 'VariableNames', {'ecRxn','stage','source','kcat','candidates','spreadLog10'});
if p.writeReport
    outDir = fullfile(params.path, 'output');
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    id = 'ecModel';
    if isfield(model, 'id') && ~isempty(model.id)
        id = char(model.id);
    end
    writetable(decisions, fullfile(outDir, [id '_kcatDecisions.tsv']), 'FileType', 'text', ...
        'Delimiter', '\t');
end
nOpen = sum(strcmp(stage, ''));
fprintf('Kcats decided for %d of %d reactions.\n', numel(stage) - nOpen, numel(stage));
end

function s = normalizeSource(label)
% As mergeKcats folds source labels: lower case, other characters as '_'.
s = lower(strtrim(char(label)));
s = regexprep(s, '[^0-9a-z]+', '_');
s = regexprep(s, '^_+|_+$', '');
end

function name = sourceFile(token)
% The source file a stage token reads: the BRENDA tiers read brenda.tsv.
if startsWith(token, 'database_')
    name = 'brenda';
else
    name = token;
end
end
