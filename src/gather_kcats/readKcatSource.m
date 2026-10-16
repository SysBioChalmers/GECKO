function [records, meta] = readKcatSource(filename)
% readKcatSource  Read a kcat source file as it is, without a model.
%
% Reads a file written by saveKcatSource or a kcatsFrom<Source> function,
% or curated by hand in the same format. Use loadKcatSource to map the
% rows to an ecModel.
%
% The file is read as UTF-8 (a leading byte-order mark is ignored). The
% '#' lines before the column header are the header (keys in lower case);
% '#' lines and blank lines further down are skipped. Fields wrapped in double quotes, as a
% spreadsheet may save them, are unwrapped. kcat and log10_sd must be plain
% numbers (12, 0.5, 1e-05); other text, such as '5,5' with a decimal comma,
% is read as NaN, which loadKcatSource reports. A row with more filled
% fields than there are columns is an error.
%
% Parameters
% ----------
% filename : char
%     path to the kcat source file.
%
% Returns
% -------
% records : struct
%     one field per column, each with one entry per data row: kcat and
%     log10_sd as double (NaN when empty or not a number), every other
%     column as a cell of char. Columns other than the standard ones are
%     kept. records.line is the line number of each row in the file.
% meta : struct
%     the '#' header lines, as field: value (for instance meta.source,
%     meta.kind).
%
% Examples
% --------
%     [records, meta] = readKcatSource(fullfile(projectPath,'data','kcats','brenda.tsv'));
%
% See also
% --------
% loadKcatSource, saveKcatSource

text = readUtf8(filename);
lines = regexp(text, '\r?\n', 'split');

meta = struct();
k = 1;
while k <= numel(lines) && (isempty(strtrim(lines{k})) || startsWith(unquoteField(lines{k}), '#'))
    [key, value] = parseHeaderLine(lines{k});
    if ~isempty(key)
        meta.(key) = value;
    end
    k = k + 1;
end
if k > numel(lines)
    error('readKcatSource:noHeader', '%s has no column header.', filename)
end
columns = cellfun(@unquoteField, strsplit(lines{k}, sprintf('\t'), ...
    'CollapseDelimiters', false), 'UniformOutput', false);
if ~all(ismember({'kcat','source'}, columns)) || ~any(ismember({'rxn','enzymes'}, columns))
    error('readKcatSource:columns', ['%s needs the columns kcat and source, and rxn ' ...
        'and/or enzymes.'], filename)
end
named = columns(~cellfun(@isempty, columns));
if numel(unique(named)) < numel(named) || any(strcmp(columns, 'line'))
    error('readKcatSource:columns', '%s has a column name twice, or a column named line.', filename)
end
nCol = numel(columns);

trimmed = strtrim(lines);
skip = cellfun(@isempty, trimmed) | startsWith(trimmed, '#') | startsWith(trimmed, '"#');
skip(1:k) = true;
kept = find(~skip);
L = string(lines(kept));
L = L(:);
F = strings(numel(kept), nCol);
nTabs = count(L, sprintf('\t'));
regular = nTabs == nCol - 1;
if any(regular)
    % Rows with exactly one field per column, split in one go.
    F(regular, :) = reshape(split(L(regular), sprintf('\t'), 2), [], nCol);
end
for r = find(~regular)'
    fields = strsplit(char(L(r)), sprintf('\t'), 'CollapseDelimiters', false);
    if numel(fields) > nCol
        if any(~cellfun(@isempty, strtrim(fields(nCol+1:end))))
            error('readKcatSource:fields', '%s, line %d has more fields than the %d columns.', ...
                filename, kept(r), nCol)
        end
        fields = fields(1:nCol);
    end
    fields(end+1:nCol) = {''};
    F(r, :) = string(fields);
end
F = strtrim(F);
quoted = find(startsWith(F, '"'));
for q = quoted(:)'
    F(q) = unquoteField(char(F(q)));
end
filled = any(F ~= "", 2); % a row of tabs only is skipped
values = cellstr(F(filled, :));
if isempty(values)
    values = cell(0, nCol);
end
lineNo = kept(filled);
lineNo = lineNo(:);
n = numel(lineNo);

records = struct();
for c = 1:nCol
    if isempty(columns{c})
        continue
    end
    name = matlab.lang.makeValidName(columns{c});
    v = values(:, c);
    if any(strcmp(columns{c}, {'kcat','log10_sd'}))
        v = parseKcatNumber(v);
    end
    records.(name) = v;
end
for c = {'rxn','enzymes','substrate','source','detail','ec_code','note'}
    if ~isfield(records, c{1})
        records.(c{1}) = repmat({''}, n, 1);
    end
end
if ~isfield(records, 'log10_sd')
    records.log10_sd = nan(n, 1);
end
records.line = lineNo;
end

function text = readUtf8(filename)
fid = fopen(filename, 'r', 'n', 'UTF-8');
if fid < 0
    error('readKcatSource:read', 'Cannot read %s', filename)
end
cleanup = onCleanup(@() fclose(fid));
text = fread(fid, [1 Inf], '*char');
text = regexprep(text, ['^' char(65279)], '');
end
