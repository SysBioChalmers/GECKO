function data = readAdapterToml(filename)
% readAdapterToml  Read a model adapter TOML file into a struct.
%
% Reads the subset of TOML that model adapter files use: comments, bare
% keys, basic and literal strings, integers, floats, booleans, arrays of
% these values (also over several lines), tables ([kegg], [a.b.c]) and
% arrays of tables ([[a]]). Anything else, such as dates, inline tables,
% multi-line strings or quoted keys, stops with an error that names the
% line.
%
% Parameters
% ----------
% filename : char
%     path to the TOML file.
%
% Returns
% -------
% data : struct
%     one field per key, with tables as nested structs and arrays of
%     tables as struct arrays. Strings are char, numbers double, booleans
%     logical, arrays of strings a 1xn cell and arrays of numbers or
%     booleans a 1xn vector. An empty array is an empty cell.
%
% Examples
% --------
%     data = readAdapterToml(fullfile(projectPath,'model_adapter.toml'));
%
% See also
% --------
% writeAdapterToml, loadAdapterParams

text = fileread(filename);
text = strrep(text, sprintf('\r\n'), newline);
lines = strsplit(text, newline, 'CollapseDelimiters', false);

data = struct();
tablePath = {};   % keys of the current table, each element {name, arrayIndex}
arrayNames = {};  % dotted names declared as arrays of tables
i = 0;
while i < numel(lines)
    i = i + 1;
    lineNo = i;
    line = stripComment(lines{i}, filename, lineNo);
    if isempty(strtrim(line))
        continue
    end
    trimmed = strtrim(line);
    if startsWith(trimmed, '[[')
        if ~endsWith(trimmed, ']]')
            tomlError(filename, lineNo, 'array of tables header must end with ]]');
        end
        keys = splitKey(trimmed(3:end-2), filename, lineNo);
        [data, tablePath] = appendTableArray(data, keys, arrayNames, filename, lineNo);
        arrayNames = union(arrayNames, {strjoin(keys, '.')});
    elseif startsWith(trimmed, '[')
        if ~endsWith(trimmed, ']')
            tomlError(filename, lineNo, 'table header must end with ]');
        end
        keys = splitKey(trimmed(2:end-1), filename, lineNo);
        tablePath = resolvePath(data, keys, arrayNames);
        data = ensureTable(data, tablePath, filename, lineNo);
    else
        eq = findUnquoted(line, '=');
        if isempty(eq)
            tomlError(filename, lineNo, 'expected key = value');
        end
        keys = splitKey(line(1:eq-1), filename, lineNo);
        valueText = strtrim(line(eq+1:end));
        % Arrays may span several lines: read on until the brackets close.
        while startsWith(valueText, '[') && bracketDepth(valueText) > 0
            if i >= numel(lines)
                tomlError(filename, lineNo, 'array is not closed');
            end
            i = i + 1;
            valueText = [valueText ' ' strtrim(stripComment(lines{i}, filename, i))]; %#ok<AGROW>
        end
        value = parseValue(valueText, filename, lineNo);
        data = setKey(data, [tablePath, keyPath(keys)], value, filename, lineNo);
    end
end
end

function line = stripComment(line, filename, lineNo)
% Removes a # comment that is not inside a string.
inBasic = false; inLiteral = false; k = 1;
while k <= numel(line)
    c = line(k);
    if inBasic
        if c == '\'
            k = k + 1;
        elseif c == '"'
            inBasic = false;
        end
    elseif inLiteral
        if c == ''''
            inLiteral = false;
        end
    elseif c == '"'
        if k + 2 <= numel(line) && strcmp(line(k:k+2), '"""')
            tomlError(filename, lineNo, 'multi-line strings are not supported');
        end
        inBasic = true;
    elseif c == ''''
        if k + 2 <= numel(line) && strcmp(line(k:k+2), '''''''')
            tomlError(filename, lineNo, 'multi-line strings are not supported');
        end
        inLiteral = true;
    elseif c == '#'
        line = line(1:k-1);
        return
    end
    k = k + 1;
end
end

function idx = findUnquoted(line, ch)
% Position of the first ch outside a string, or empty.
idx = [];
inBasic = false; inLiteral = false; k = 1;
while k <= numel(line)
    c = line(k);
    if inBasic
        if c == '\'
            k = k + 1;
        elseif c == '"'
            inBasic = false;
        end
    elseif inLiteral
        if c == ''''
            inLiteral = false;
        end
    elseif c == '"'
        inBasic = true;
    elseif c == ''''
        inLiteral = true;
    elseif c == ch
        idx = k;
        return
    end
    k = k + 1;
end
end

function depth = bracketDepth(text)
% Open [ minus closed ] outside strings.
depth = 0; inBasic = false; inLiteral = false; k = 1;
while k <= numel(text)
    c = text(k);
    if inBasic
        if c == '\'
            k = k + 1;
        elseif c == '"'
            inBasic = false;
        end
    elseif inLiteral
        if c == ''''
            inLiteral = false;
        end
    elseif c == '"'
        inBasic = true;
    elseif c == ''''
        inLiteral = true;
    elseif c == '['
        depth = depth + 1;
    elseif c == ']'
        depth = depth - 1;
    end
    k = k + 1;
end
end

function keys = splitKey(text, filename, lineNo)
% Dotted bare key into its parts.
keys = strtrim(strsplit(strtrim(text), '.'));
for k = 1:numel(keys)
    if isempty(regexp(keys{k}, '^[A-Za-z][A-Za-z0-9_]*$', 'once'))
        tomlError(filename, lineNo, sprintf(['key "%s" is not supported: use letters, ' ...
            'digits and underscores, starting with a letter'], keys{k}));
    end
end
end

function p = keyPath(keys)
p = cellfun(@(k) {k, 0}, keys, 'UniformOutput', false);
end

function path = resolvePath(data, keys, arrayNames)
% Path of {name, index} pairs; a key that names an array of tables refers
% to its last element.
path = cell(1, numel(keys));
for k = 1:numel(keys)
    idx = 0;
    if ismember(strjoin(keys(1:k), '.'), arrayNames)
        idx = numel(getKey(data, [path(1:k-1), {{keys{k}, 0}}]));
    end
    path{k} = {keys{k}, idx};
end
end

function [data, tablePath] = appendTableArray(data, keys, arrayNames, filename, lineNo)
% Appends a new element to the array of tables named by keys.
parent = resolvePath(data, keys(1:end-1), arrayNames);
data = ensureTable(data, parent, filename, lineNo);
name = strjoin(keys, '.');
existing = getKey(data, [parent, {{keys{end}, 0}}]);
if isempty(existing) || ~ismember(name, arrayNames)
    if ~isempty(existing)
        tomlError(filename, lineNo, sprintf('"%s" is already a table or value', name));
    end
    n = 1;
else
    n = numel(existing) + 1;
end
tablePath = [parent, {{keys{end}, n}}];
data = setTableElement(data, parent, keys{end}, n);
end

function data = ensureTable(data, path, filename, lineNo)
if isempty(path)
    return
end
name = path{1}{1}; idx = path{1}{2};
if ~isfield(data, name)
    data.(name) = struct();
elseif ~isstruct(data.(name))
    tomlError(filename, lineNo, sprintf('"%s" is already a value, not a table', name));
end
child = data.(name);
if idx == 0
    child = ensureTable(child, path(2:end), filename, lineNo);
    data.(name) = child;
else
    element = child(idx);
    element = ensureTable(element, path(2:end), filename, lineNo);
    data.(name) = assignElement(child, idx, element);
end
end

function data = setTableElement(data, parent, name, n)
% Adds element n (an empty struct) to the table array name under parent.
if isempty(parent)
    if ~isfield(data, name) || isempty(data.(name)) || n == 1
        data.(name) = struct();
    else
        data.(name) = assignElement(data.(name), n, struct());
    end
    return
end
pname = parent{1}{1}; pidx = parent{1}{2};
child = data.(pname);
if pidx == 0
    data.(pname) = setTableElement(child, parent(2:end), name, n);
else
    element = setTableElement(child(pidx), parent(2:end), name, n);
    data.(pname) = assignElement(child, pidx, element);
end
end

function value = getKey(data, path)
value = [];
node = data;
for k = 1:numel(path)
    name = path{k}{1}; idx = path{k}{2};
    if ~isstruct(node) || ~isfield(node, name)
        value = [];
        return
    end
    node = node.(name);
    if idx > 0
        node = node(idx);
    end
end
value = node;
end

function data = setKey(data, path, value, filename, lineNo)
name = path{1}{1}; idx = path{1}{2};
if numel(path) == 1
    % [] only fills a field that another element of a table array has;
    % TOML itself has no empty value.
    if isfield(data, name) && ~isequal(data.(name), [])
        tomlError(filename, lineNo, sprintf('key "%s" is defined twice', name));
    end
    data.(name) = value;
    return
end
if ~isfield(data, name)
    data.(name) = struct();
elseif ~isstruct(data.(name))
    tomlError(filename, lineNo, sprintf('"%s" is already a value, not a table', name));
end
child = data.(name);
if idx == 0
    data.(name) = setKey(child, path(2:end), value, filename, lineNo);
else
    element = setKey(child(idx), path(2:end), value, filename, lineNo);
    data.(name) = assignElement(child, idx, element);
end
end

function s = assignElement(s, idx, element)
% Assigns element to s(idx), adding fields so the struct array stays
% consistent.
if idx > numel(s) || isempty(fieldnames(s))
    if numel(s) == 1 && isempty(fieldnames(s)) && idx == 1
        s = element;
        return
    end
end
fs = fieldnames(element);
for f = 1:numel(fs)
    if ~isfield(s, fs{f})
        [s.(fs{f})] = deal([]);
    end
end
fs = fieldnames(s);
for f = 1:numel(fs)
    if ~isfield(element, fs{f})
        element.(fs{f}) = [];
    end
end
element = orderfields(element, s(1));
s(idx) = element;
end

function value = parseValue(text, filename, lineNo)
text = strtrim(text);
if isempty(text)
    tomlError(filename, lineNo, 'missing value');
end
switch text(1)
    case '"'
        [value, rest] = parseBasicString(text, filename, lineNo);
    case ''''
        close = find(text(2:end) == '''', 1) + 1;
        if isempty(close)
            tomlError(filename, lineNo, 'string is not closed');
        end
        value = text(2:close-1);
        rest = text(close+1:end);
    case '['
        [value, rest] = parseArray(text, filename, lineNo);
    case '{'
        tomlError(filename, lineNo, 'inline tables are not supported; use a [table] instead');
    otherwise
        [value, rest] = parseScalar(text, filename, lineNo);
end
if ~isempty(strtrim(rest))
    tomlError(filename, lineNo, sprintf('unexpected text after the value: %s', strtrim(rest)));
end
end

function [value, rest] = parseBasicString(text, filename, lineNo)
value = '';
k = 2;
while k <= numel(text)
    c = text(k);
    if c == '"'
        rest = text(k+1:end);
        return
    elseif c == '\'
        if k == numel(text)
            break
        end
        e = text(k+1);
        switch e
            case 'n', value(end+1) = newline;
            case 't', value(end+1) = sprintf('\t');
            case '"', value(end+1) = '"';
            case '\', value(end+1) = '\';
            case 'u'
                value(end+1) = char(hex2dec(text(k+2:k+5)));
                k = k + 4;
            otherwise
                tomlError(filename, lineNo, sprintf('escape \\%s is not supported', e));
        end
        k = k + 2;
        continue
    else
        value(end+1) = c; %#ok<AGROW>
    end
    k = k + 1;
end
tomlError(filename, lineNo, 'string is not closed');
end

function [value, rest] = parseArray(text, filename, lineNo)
% Parses [a, b, ...]; returns a cell for strings, a vector otherwise.
items = {};
k = 2;
while true
    while k <= numel(text) && (isspace(text(k)) || text(k) == ',')
        k = k + 1;
    end
    if k > numel(text)
        tomlError(filename, lineNo, 'array is not closed');
    end
    if text(k) == ']'
        rest = text(k+1:end);
        break
    end
    if text(k) == '['
        tomlError(filename, lineNo, 'nested arrays are not supported');
    end
    % Find the end of this element.
    if text(k) == '"'
        [item, after] = parseBasicString(text(k:end), filename, lineNo);
        k = numel(text) - numel(after) + 1;
    elseif text(k) == ''''
        close = find(text(k+1:end) == '''', 1) + k;
        if isempty(close)
            tomlError(filename, lineNo, 'string is not closed');
        end
        item = text(k+1:close-1);
        k = close + 1;
    else
        stop = k;
        while stop <= numel(text) && text(stop) ~= ',' && text(stop) ~= ']'
            stop = stop + 1;
        end
        [item, after] = parseScalar(strtrim(text(k:stop-1)), filename, lineNo);
        if ~isempty(strtrim(after))
            tomlError(filename, lineNo, 'unexpected text in array');
        end
        k = stop;
    end
    items{end+1} = item; %#ok<AGROW>
end
if isempty(items)
    value = {};
elseif all(cellfun(@ischar, items))
    value = items;
elseif all(cellfun(@(x) isnumeric(x) && isscalar(x), items))
    value = cell2mat(items);
elseif all(cellfun(@islogical, items))
    value = cell2mat(items);
else
    tomlError(filename, lineNo, 'arrays must hold values of one type');
end
end

function [value, rest] = parseScalar(text, filename, lineNo)
token = regexp(text, '^[^\s,\]]+', 'match', 'once');
rest = text(numel(token)+1:end);
switch token
    case 'true'
        value = true;
    case 'false'
        value = false;
    case {'inf', '+inf'}
        value = Inf;
    case '-inf'
        value = -Inf;
    case {'nan', '+nan', '-nan'}
        value = NaN;
    otherwise
        clean = strrep(token, '_', '');
        if ~isempty(regexp(clean, '^[+-]?(0|[1-9][0-9]*)$', 'once')) || ...
                ~isempty(regexp(clean, '^[+-]?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?$', 'once'))
            value = str2double(clean);
        elseif ~isempty(regexp(token, '^\d{4}-\d{2}-\d{2}', 'once'))
            tomlError(filename, lineNo, 'dates are not supported');
        else
            tomlError(filename, lineNo, sprintf('cannot read value "%s"; strings need quotes', token));
        end
end
end

function tomlError(filename, lineNo, message)
error('GECKO:readAdapterToml', '%s, line %d: %s', filename, lineNo, message);
end
