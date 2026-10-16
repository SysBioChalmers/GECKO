function writeAdapterToml(data, filename, header)
% writeAdapterToml  Write a struct as a model adapter TOML file.
%
% Writes the subset of TOML that readAdapterToml reads: values first,
% then each nested struct as a [table] and each struct array as an array
% of tables ([[name]]).
%
% Parameters
% ----------
% data : struct
%     values may be char or string scalars, numeric or logical scalars or
%     vectors, cell arrays of char, structs and struct arrays.
% filename : char
%     path of the TOML file to write; an existing file is overwritten.
% header : char or cell, optional
%     comment lines written at the top of the file (without '#').
%
% Examples
% --------
%     writeAdapterToml(struct('sigma',0.5), 'model_adapter.toml', ...
%         {'Adapter for ecYeastGEM'});
%
% See also
% --------
% readAdapterToml, convertAdapterToToml

if nargin < 3
    header = {};
elseif ischar(header) || isstring(header)
    header = cellstr(header);
end

out = {};
for k = 1:numel(header)
    out{end+1} = ['# ' header{k}]; %#ok<AGROW>
end
if ~isempty(header)
    out{end+1} = '';
end
out = writeTable(out, data, {});

fid = fopen(filename, 'w');
if fid < 0
    error('GECKO:writeAdapterToml', 'Cannot write %s', filename);
end
fprintf(fid, '%s\n', out{:});
fclose(fid);
end

function out = writeTable(out, s, prefix)
fields = fieldnames(s);
tables = {};
for k = 1:numel(fields)
    value = s.(fields{k});
    if isstruct(value)
        tables{end+1} = fields{k}; %#ok<AGROW>
    elseif isempty(value) && ~iscell(value) && ~ischar(value)
        continue % unset value; TOML has no null
    else
        out{end+1} = [fields{k} ' = ' formatValue(value, strjoin([prefix fields(k)], '.'))]; %#ok<AGROW>
    end
end
for k = 1:numel(tables)
    value = s.(tables{k});
    name = strjoin([prefix tables(k)], '.');
    if numel(value) == 1
        out{end+1} = ''; %#ok<AGROW>
        out{end+1} = ['[' name ']']; %#ok<AGROW>
        out = writeTable(out, value, [prefix tables(k)]);
    else
        for e = 1:numel(value)
            out{end+1} = ''; %#ok<AGROW>
            out{end+1} = ['[[' name ']]']; %#ok<AGROW>
            out = writeTable(out, value(e), [prefix tables(k)]);
        end
    end
end
end

function text = formatValue(value, name)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ischar(value)
    text = quote(value);
elseif iscell(value)
    if ~all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), value))
        error('GECKO:writeAdapterToml', '%s: only cell arrays of text can be written', name);
    end
    parts = cellfun(@(x) quote(char(x)), value, 'UniformOutput', false);
    text = ['[' strjoin(parts(:)', ', ') ']'];
elseif islogical(value) || isnumeric(value)
    if isscalar(value)
        text = formatScalar(value);
    else
        parts = arrayfun(@formatScalar, value(:)', 'UniformOutput', false);
        text = ['[' strjoin(parts, ', ') ']'];
    end
else
    error('GECKO:writeAdapterToml', '%s: values of class %s cannot be written', name, class(value));
end
end

function text = formatScalar(x)
if islogical(x)
    if x
        text = 'true';
    else
        text = 'false';
    end
elseif isinf(x)
    if x > 0
        text = 'inf';
    else
        text = '-inf';
    end
elseif isnan(x)
    text = 'nan';
elseif x == round(x) && abs(x) < 1e15
    text = sprintf('%d', x);
else
    text = sprintf('%.17g', x);
    % Use the shortest form that reads back to the same double.
    for p = 1:17
        candidate = sprintf('%.*g', p, x);
        if str2double(candidate) == x
            text = candidate;
            break
        end
    end
    if isempty(regexp(text, '[.eE]', 'once'))
        text = [text '.0'];
    end
end
end

function text = quote(s)
s = strrep(s, '\', '\\');
s = strrep(s, '"', '\"');
s = strrep(s, sprintf('\t'), '\t');
s = strrep(s, newline, '\n');
text = ['"' s '"'];
end
