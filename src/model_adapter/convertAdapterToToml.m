function report = convertAdapterToToml(adapter, tomlFile)
% convertAdapterToToml  Write the parameters of an adapter class to model_adapter.toml.
%
% Converts an existing ModelAdapter class (GECKO 3 or 4 style) to
% model_adapter.toml, the configuration file shared by GECKO and geckopy.
% Parameter names are translated with adapterParameters.tsv (for instance
% Ptot becomes p_tot), paths inside the project folder are written
% relative to it, and parameters without a geckopy name are kept under
% [matlab.params], so that loadAdapterParams gives back the same params.
%
% Methods that the adapter class overrides cannot be written to TOML.
% They are listed in the report, with how to carry them over:
% - getSpontaneousReactions: set [spontaneous] name_contains or file;
% - makeModelAnaerobic, changeProteinBiomass: set [matlab]
%   make_model_anaerobic or change_protein_biomass to a function name;
% - other methods: keep the adapter class and let its constructor call
%   loadAdapterParams.
%
% Parameters
% ----------
% adapter : ModelAdapter or char
%     an adapter object, or the path to its .m file.
% tomlFile : char, optional
%     file to write (default: model_adapter.toml in the adapter's
%     params.path). An existing file is not overwritten.
%
% Returns
% -------
% report : struct
%     - file: the file written;
%     - unmapped: parameters written under [matlab.params];
%     - methods: overridden methods that need carrying over by hand.
%
% Examples
% --------
%     report = convertAdapterToToml(fullfile(projectPath,'MyGEMAdapter.m'));
%     adapter = TomlAdapter(projectPath);
%
% See also
% --------
% TomlAdapter, loadAdapterParams, writeAdapterToml

if ischar(adapter) || isstring(adapter)
    adapter = ModelAdapterManager.getAdapter(char(adapter));
end
params = adapter.params;
if isfield(params,'path') && isstruct(params.path)
    params.path = params.path.Name; % a struct path, as some old adapters keep it
end
projectPath = '';
if isfield(params,'path')
    projectPath = char(params.path);
end
if nargin < 2 || isempty(tomlFile)
    tomlFile = fullfile(projectPath, 'model_adapter.toml');
end
if isfile(tomlFile)
    error('GECKO:convertAdapterToToml', '%s exists already; it is not overwritten.', tomlFile);
end

map = readMap();
data = struct();
used = {};
for k = 1:numel(map.tomlKey)
    field = map.matlabField{k};
    if isempty(field) || any(strcmp(field, {'evotune.kcatSources','evotune.sigma0logSource'}))
        continue
    end
    [found, value] = getDotted(params, field);
    if ~found || (isempty(value) && ~ischar(value))
        if found
            used{end+1} = field; %#ok<AGROW>
        end
        continue
    end
    used{end+1} = field; %#ok<AGROW>
    if any(strcmp(map.tomlKey{k}, {'conv_gem','spontaneous.file'}))
        value = relativePath(char(value), projectPath);
    end
    data = setDotted(data, map.tomlKey{k}, value);
end
used{end+1} = 'path';

% Trust tiers: one group per ec.source value.
if isfield(params,'evotune') && isfield(params.evotune,'kcatSources')
    used = [used, {'evotune.kcatSources','evotune.sigma0logSource'}];
    sources = params.evotune.kcatSources;
    widths = params.evotune.sigma0logSource;
    groups = struct();
    sigma0 = struct();
    for k = 1:numel(sources)
        name = matlab.lang.makeValidName(sources{k});
        groups.(name).sources = sources(k);
        sigma0.(name) = widths(k);
    end
    data.evotune.source_groups = groups;
    data.evotune.sigma0_log_source = sigma0;
end

% Everything else is kept for GECKO under [matlab.params].
extra = unusedFields(params, '', used);
report.unmapped = {};
if ~isempty(fieldnames(extra))
    data.matlab.params = extra;
    report.unmapped = flattenNames(extra, '');
end

report.methods = overriddenMethods(adapter);
report.file = tomlFile;

header = {sprintf('Model adapter for %s, converted from %s.', ...
    getOr(params,'org_name','this project'), class(adapter)), ...
    'Shared by GECKO (MATLAB) and geckopy; see adapterParameters.tsv in GECKO for every key.'};
writeAdapterToml(data, tomlFile, header);

fprintf('Wrote %s.\n', tomlFile);
if ~isempty(report.unmapped)
    fprintf('Kept for GECKO only, under [matlab.params]: %s\n', strjoin(report.unmapped, ', '));
end
if ~isempty(report.methods)
    printOrange(sprintf(['%s overrides methods that a TOML file cannot hold: %s.\n' ...
        'Set [spontaneous] for getSpontaneousReactions, [matlab] function names for\n' ...
        'makeModelAnaerobic and changeProteinBiomass, or keep the class and let its\n' ...
        'constructor call loadAdapterParams.\n'], class(adapter), strjoin(report.methods, ', ')));
end
end

function map = readMap()
file = fullfile(fileparts(mfilename('fullpath')), 'adapterParameters.tsv');
fid = fopen(file, 'r');
cols = textscan(fid, '%s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1, 'Whitespace', '');
fclose(fid);
map.tomlKey = cols{1};
map.matlabField = cols{2};
end

function names = overriddenMethods(adapter)
% Methods defined by the adapter class itself rather than ModelAdapter.
names = {};
mc = metaclass(adapter);
for k = 1:numel(mc.MethodList)
    m = mc.MethodList(k);
    if strcmp(m.DefiningClass.Name, mc.Name) && ~strcmp(m.Name, mc.Name) && ~m.Static
        names{end+1} = m.Name; %#ok<AGROW>
    end
end
names = sort(names);
end

function extra = unusedFields(s, prefix, used)
extra = struct();
fields = fieldnames(s);
for k = 1:numel(fields)
    name = fields{k};
    if ~isempty(prefix)
        name = [prefix '.' name];
    end
    if any(strcmp(name, used))
        continue
    end
    value = s.(fields{k});
    if isstruct(value) && isscalar(value) && any(startsWith(used, [name '.']))
        sub = unusedFields(value, name, used);
        if ~isempty(fieldnames(sub))
            extra.(fields{k}) = sub;
        end
    elseif ~(isempty(value) && ~ischar(value) && ~iscell(value))
        extra.(fields{k}) = value;
    end
end
end

function names = flattenNames(s, prefix)
names = {};
fields = fieldnames(s);
for k = 1:numel(fields)
    name = fields{k};
    if ~isempty(prefix)
        name = [prefix '.' name];
    end
    if isstruct(s.(fields{k}))
        names = [names, flattenNames(s.(fields{k}), name)]; %#ok<AGROW>
    else
        names{end+1} = name; %#ok<AGROW>
    end
end
end

function p = relativePath(p, projectPath)
p = strrep(p, '\', '/');
base = strrep(projectPath, '\', '/');
if ~isempty(base) && startsWith(p, [base '/'])
    p = p(numel(base)+2:end);
end
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

function v = getOr(s, field, default)
v = default;
if isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
end
end
