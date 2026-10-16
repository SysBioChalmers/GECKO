function config = kcatAssignmentConfig(params)
% kcatAssignmentConfig  The kcat assignment settings of a model adapter.
%
% Reads [kcat_assignment] and [kcat_sources.*] of model_adapter.toml (or
% params.kcat_assignment and params.kcat_sources of an adapter class), and
% expands a preset into its stages. Without [kcat_assignment] the preset
% is gecko_classic, which reproduces the GECKO tutorial.
%
% Returns a struct with:
% - preset: the preset name, or '' when stages are given;
% - isozymes, standard: logical;
% - stages: struct array with action, sources (cell), criteria,
%   topOriginLimit, bottomOriginLimit, wildcardLimit;
% - sources: struct, per source name its options (from [kcat_sources]).
% geckopy's kcat_assignment_config gives the same settings.

a = struct();
if isfield(params, 'kcat_assignment') && isstruct(params.kcat_assignment)
    a = params.kcat_assignment;
end
config.preset = '';
config.isozymes = getOr(a, 'isozymes', true);
config.standard = getOr(a, 'standard', true);
config.sources = struct();
hasPreset = isfield(a, 'preset') && ~isempty(a.preset);
hasStages = isfield(a, 'stage') && ~isempty(a.stage);
if hasPreset && hasStages
    error('kcatAssignment:config', ['[kcat_assignment] gives both a preset and stages; ' ...
        'give one of the two.'])
end
if hasStages
    config.stages = arrayfun(@readStage, a.stage(:)');
else
    preset = 'gecko_classic';
    if hasPreset
        preset = char(a.preset);
    end
    [config.stages, config.sources] = presetStages(preset);
    config.preset = preset;
end
if isfield(params, 'kcat_sources') && isstruct(params.kcat_sources)
    names = fieldnames(params.kcat_sources);
    for k = 1:numel(names)
        config.sources.(names{k}) = params.kcat_sources.(names{k});
    end
end
end

function [stages, sources] = presetStages(preset)
% The presets of kcat assignment; geckopy defines the same.
select = @(srcs) struct('action', 'select', 'sources', {srcs}, 'criteria', 'max', ...
    'topOriginLimit', 6, 'bottomOriginLimit', 6, 'wildcardLimit', 3);
switch preset
    case 'gecko_classic'
        stages = select({'custom', 'database_top', 'dlkcat', 'database_bottom'});
        sources.brenda = struct('aggregation', 'max');
    case 'gecko_median'
        stages = select({'custom', 'database_top', 'dlkcat', 'database_bottom'});
        sources.brenda = struct('aggregation', 'median');
    otherwise
        error('kcatAssignment:config', 'Unknown kcat assignment preset "%s". Presets: %s.', ...
            preset, 'gecko_classic, gecko_median')
end
end

function stage = readStage(s)
stage.action = getOr(s, 'action', 'select');
if ~strcmp(stage.action, 'select')
    error('kcatAssignment:config', ['Unknown stage action "%s"; this GECKO version knows ' ...
        '"select".'], stage.action)
end
sources = getOr(s, 'sources', {});
if ischar(sources)
    sources = {sources};
end
if isempty(sources)
    error('kcatAssignment:config', 'A [[kcat_assignment.stage]] needs sources.')
end
stage.sources = sources(:)';
stage.criteria = validatestring(getOr(s, 'criteria', 'max'), {'max','min','median','mean'});
stage.topOriginLimit = getOr(s, 'top_origin_limit', 6);
stage.bottomOriginLimit = getOr(s, 'bottom_origin_limit', 6);
stage.wildcardLimit = getOr(s, 'wildcard_limit', 3);
end

function v = getOr(s, field, default)
v = default;
if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
end
end
