function params = loadAdapterParams(tomlFile)
% loadAdapterParams  Read model adapter parameters from a TOML file.
%
% Reads model_adapter.toml, the configuration file shared by GECKO and
% geckopy, and returns the parameters with GECKO's field names, so that
% they can be used as ModelAdapter.params. Keys are translated with
% adapterParameters.tsv (for instance p_tot becomes Ptot); keys that are
% not given get the default listed there, which is geckopy's default.
%
% Parameters
% ----------
% tomlFile : char
%     path to model_adapter.toml, or to the project folder that holds it.
%
% Returns
% -------
% params : struct
%     the adapter parameters, as in ModelAdapter.params:
%     - path: the project folder (the folder of the TOML file);
%     - convGEM and spontaneous.file: absolute paths;
%     - evotune.kcatSources and evotune.sigma0logSource: built from the
%       evotune.source_groups and evotune.sigma0_log_source tables, one
%       entry per ec.source value;
%     - fields under [matlab.params] are copied as they are, for
%       parameters that only GECKO uses.
%
% Examples
% --------
%     params = loadAdapterParams(fullfile(projectPath,'model_adapter.toml'));
%
% See also
% --------
% TomlAdapter, readAdapterToml, convertAdapterToToml

if isfolder(tomlFile)
    tomlFile = fullfile(tomlFile, 'model_adapter.toml');
end
if ~isfile(tomlFile)
    error('GECKO:loadAdapterParams', 'No model adapter file found at %s', tomlFile);
end
tomlFile = char(tomlFile);
[projectPath, ~, ~] = fileparts(getAbsolutePath(tomlFile));

data = readAdapterToml(tomlFile);
map = adapterParameterMap();

% Unknown keys are errors, as in geckopy, so that a typo is not ignored.
known = [map.tomlKey; {'evotune.source_groups'; 'evotune.sigma0_log_source'}];
checkKnownKeys(data, '', known, tomlFile);

params = struct();
params.path = projectPath;
for k = 1:numel(map.tomlKey)
    key = map.tomlKey{k};
    if any(strcmp(key, {'evotune.source_groups', 'evotune.sigma0_log_source'}))
        continue
    end
    [found, value] = getDotted(data, key);
    if ~found
        if isempty(map.default{k})
            continue
        end
        value = map.default{k};
    end
    if isempty(map.matlabField{k})
        field = key; % a key that only geckopy reads keeps its name
    else
        field = map.matlabField{k};
    end
    params = setDotted(params, field, value);
end

if ~isfield(params, 'complex')
    params.complex.taxonomicID = []; % geckopy's default: no Complex Portal query
end

for f = {'convGEM'}
    if isfield(params, f{1})
        params.(f{1}) = projectFile(projectPath, params.(f{1}));
    end
end
if isfield(params, 'spontaneous') && isfield(params.spontaneous, 'file')
    params.spontaneous.file = projectFile(projectPath, params.spontaneous.file);
end

params = sourceGroups(params, data, tomlFile);

if isfield(data, 'matlab') && isfield(data.matlab, 'params')
    params = mergeStruct(params, data.matlab.params);
end
end

function map = adapterParameterMap()
% Reads adapterParameters.tsv into toml keys, MATLAB fields and defaults.
file = fullfile(fileparts(mfilename('fullpath')), 'adapterParameters.tsv');
fid = fopen(file, 'r');
cols = textscan(fid, '%s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1, ...
    'Whitespace', '', 'ReturnOnError', false);
fclose(fid);
map.tomlKey = cols{1};
map.matlabField = cols{2};
map.default = cell(size(cols{3}));
for k = 1:numel(cols{3})
    if ~isempty(cols{3}{k})
        map.default{k} = parseDefault(cols{3}{k});
    end
end
end

function value = parseDefault(text)
% A default is written as a TOML value.
tmp = [tempname '.toml'];
fid = fopen(tmp, 'w');
fprintf(fid, 'x = %s\n', text);
fclose(fid);
cleanup = onCleanup(@() delete(tmp));
data = readAdapterToml(tmp);
value = data.x;
if ischar(value) && isempty(value)
    value = '';
end
end

function checkKnownKeys(s, prefix, known, tomlFile)
fields = fieldnames(s);
for k = 1:numel(fields)
    key = fields{k};
    if ~isempty(prefix)
        key = [prefix '.' key];
    end
    if strcmp(key, 'matlab.params')
        continue % free-form parameters for GECKO only
    end
    if any(strcmp(key, known)) || any(startsWith(known, [key '.']))
        if isstruct(s.(fields{k})) && ~any(strcmp(key, known))
            checkKnownKeys(s.(fields{k}), key, known, tomlFile);
        end
        continue
    end
    if any(startsWith(key, {'evotune.source_groups.', 'evotune.sigma0_log_source.'}))
        continue
    end
    error('GECKO:loadAdapterParams', ['Unknown key "%s" in %s. Parameters that only ' ...
        'GECKO uses go under [matlab.params].'], key, tomlFile);
end
end

function params = sourceGroups(params, data, tomlFile)
% evotune.source_groups and evotune.sigma0_log_source become the
% kcatSources and sigma0logSource lists that evotune reads: one entry per
% ec.source value, with its group's prior width.
groups = struct();
widths = struct();
if isfield(data, 'evotune') && isfield(data.evotune, 'source_groups')
    groups = data.evotune.source_groups;
    if isfield(data.evotune, 'sigma0_log_source')
        widths = data.evotune.sigma0_log_source;
    end
else
    % geckopy's defaults (EvotuneParams.source_groups, sigma0_log_source).
    groups.dlkcat.sources = {'dlkcat'};
    groups.brenda.sources = {'brenda'};
    groups.custom.sources = {'custom'};
    widths.dlkcat = 0.4;
    widths.brenda = 0.2;
    widths.custom = 0.1;
end
names = fieldnames(groups);
if ~isempty(setxor(names, fieldnames(widths)))
    error('GECKO:loadAdapterParams', ['%s: evotune.sigma0_log_source must have ' ...
        'the same group names as evotune.source_groups'], tomlFile);
end
kcatSources = {};
sigma0logSource = [];
for k = 1:numel(names)
    group = groups.(names{k});
    sources = {};
    if isfield(group, 'sources')
        sources = group.sources;
        if ischar(sources)
            sources = {sources};
        end
    end
    if isfield(group, 'match_okp') && group.match_okp && isfield(params, 'okp')
        sources{end+1} = params.okp.method; %#ok<AGROW>
    end
    kcatSources = [kcatSources, sources]; %#ok<AGROW>
    sigma0logSource = [sigma0logSource; repmat(widths.(names{k}), numel(sources), 1)]; %#ok<AGROW>
end
params.evotune.kcatSources = kcatSources;
params.evotune.sigma0logSource = sigma0logSource;
end

function [found, value] = getDotted(s, key)
parts = strsplit(key, '.');
value = s;
found = true;
for k = 1:numel(parts)
    if ~isstruct(value) || ~isfield(value, parts{k})
        found = false;
        value = [];
        return
    end
    value = value.(parts{k});
end
end

function s = setDotted(s, key, value)
parts = strsplit(key, '.');
if numel(parts) == 1
    s.(parts{1}) = value;
    return
end
if ~isfield(s, parts{1}) || ~isstruct(s.(parts{1}))
    s.(parts{1}) = struct();
end
s.(parts{1}) = setDotted(s.(parts{1}), strjoin(parts(2:end), '.'), value);
end

function s = mergeStruct(s, extra)
fields = fieldnames(extra);
for k = 1:numel(fields)
    value = extra.(fields{k});
    if isstruct(value) && isscalar(value) && isfield(s, fields{k}) && isstruct(s.(fields{k}))
        s.(fields{k}) = mergeStruct(s.(fields{k}), value);
    else
        s.(fields{k}) = value;
    end
end
end

function p = projectFile(projectPath, p)
if isempty(p)
    return
end
if ~isAbsolute(p)
    p = fullfile(projectPath, p);
end
end

function tf = isAbsolute(p)
tf = startsWith(p, '/') || startsWith(p, '\') || ~isempty(regexp(p, '^[A-Za-z]:', 'once'));
end

function p = getAbsolutePath(p)
if ~isAbsolute(p)
    p = fullfile(pwd, p);
end
end
