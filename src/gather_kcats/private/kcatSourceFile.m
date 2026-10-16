function filename = kcatSourceFile(name, folder, modelAdapter)
% kcatSourceFile  Path of a kcat source file that may be written.
%
% Used by saveKcatSource and the kcatsFrom<Source> functions before any
% work is done: checks the name, resolves the folder (data/kcats in the
% project folder by default), and refuses a file that is curated by hand:
% named custom.tsv, or whose header says 'kind: curated' or 'source:
% custom' (in any case, also when a spreadsheet quoted the line).

name = char(name);
if isempty(regexp(name, '^[A-Za-z0-9_.-]+$', 'once'))
    error('saveKcatSource:invalidName', ...
        'The source name "%s" may only contain letters, digits, ''_'', ''-'' and ''.''.', name)
end
if isempty(folder)
    if isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    folder = fullfile(modelAdapter.params.path, 'data', 'kcats');
end
filename = fullfile(char(folder), [name '.tsv']);
if strcmpi(name, 'custom')
    error('saveKcatSource:curated', ['custom.tsv holds curated kcats and is not ' ...
        'written by GECKO; edit it by hand.'])
end
if isfile(filename)
    meta = readHeaderLines(filename);
    if (isfield(meta, 'kind') && strcmpi(meta.kind, 'curated')) || ...
            (isfield(meta, 'source') && strcmpi(meta.source, 'custom'))
        error('saveKcatSource:curated', '%s is curated by hand and is not overwritten.', filename)
    end
end
end

function meta = readHeaderLines(filename)
% Only the header lines at the top, so that a damaged table does not stop
% a gathered file from being written again.
meta = struct();
fid = fopen(filename, 'r', 'n', 'UTF-8');
if fid < 0
    return
end
cleanup = onCleanup(@() fclose(fid));
line = fgetl(fid);
if ischar(line)
    line = regexprep(line, ['^' char(65279)], '');
end
while ischar(line) && (isempty(strtrim(line)) || startsWith(unquoteField(line), '#'))
    [key, value] = parseHeaderLine(line);
    if ~isempty(key)
        meta.(key) = value;
    end
    line = fgetl(fid);
end
end
