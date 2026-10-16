function j = kcatSourceJson(v)
% kcatSourceJson  The parameters line of kcat source files as JSON.
%
% Compact JSON with keys sorted at every level. Numbers are written with
% formatKcatNumber, and NaN and infinity as null; logical values as true
% or false; char as a string (non-ASCII characters as they are); a cell
% array as a list; a numeric or logical vector with more than one element
% as a list. geckopy's save_kcat_source writes the same text.
if isstruct(v) && isscalar(v)
    f = sort(fieldnames(v));
    parts = cell(1, numel(f));
    for k = 1:numel(f)
        parts{k} = [jsonString(f{k}) ':' kcatSourceJson(v.(f{k}))];
    end
    j = ['{' strjoin(parts, ',') '}'];
elseif ischar(v) || (isstring(v) && isscalar(v))
    j = jsonString(char(v));
elseif iscell(v)
    j = ['[' strjoin(cellfun(@kcatSourceJson, v(:)', 'UniformOutput', false), ',') ']'];
elseif (isnumeric(v) || islogical(v)) && numel(v) ~= 1
    j = ['[' strjoin(arrayfun(@kcatSourceJson, v(:)', 'UniformOutput', false), ',') ']'];
elseif islogical(v)
    j = 'false';
    if v
        j = 'true';
    end
elseif isnumeric(v)
    if ~isfinite(v)
        j = 'null';
    else
        j = formatKcatNumber(double(v));
    end
elseif isempty(v)
    j = 'null';
else
    error('kcatSourceJson:type', 'Parameters of class %s cannot be written.', class(v))
end
end

function s = jsonString(s)
% JSON escapes, as Python's json.dumps(ensure_ascii=False).
out = '';
for c = s
    switch c
        case '"', out = [out '\"']; %#ok<AGROW>
        case '\', out = [out '\\']; %#ok<AGROW>
        case char(8), out = [out '\b']; %#ok<AGROW>
        case char(9), out = [out '\t']; %#ok<AGROW>
        case char(10), out = [out '\n']; %#ok<AGROW>
        case char(12), out = [out '\f']; %#ok<AGROW>
        case char(13), out = [out '\r']; %#ok<AGROW>
        otherwise
            if double(c) < 32
                out = [out sprintf('\\u%04x', double(c))]; %#ok<AGROW>
            else
                out = [out c]; %#ok<AGROW>
            end
    end
end
s = ['"' out '"'];
end
