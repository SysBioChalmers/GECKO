function [key, value] = parseHeaderLine(line)
% parseHeaderLine  Key and value of a '# key: value' header line.
%
% The line may be wrapped in quotes by a spreadsheet. The key is returned
% in lower case; it must start with a letter and be at most 63 letters,
% digits or underscores long. Other lines give an empty key. geckopy reads
% header lines the same way.
key = ''; value = '';
tok = regexp(unquoteField(line), '^\s*#\s*([A-Za-z][A-Za-z0-9_]{0,62})\s*:\s?(.*)$', 'tokens', 'once');
if ~isempty(tok)
    key = lower(tok{1});
    value = strtrim(tok{2});
end
end
