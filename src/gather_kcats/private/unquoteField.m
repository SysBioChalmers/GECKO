function t = unquoteField(t)
% unquoteField  A field without surrounding spaces and double quotes.
%
% A field wrapped in double quotes, as a spreadsheet may save it, is
% unwrapped, with "" for a quote inside; spaces inside the quotes are
% trimmed too.
t = strtrim(t);
if numel(t) >= 2 && t(1) == '"' && t(end) == '"'
    t = strtrim(strrep(t(2:end-1), '""', '"'));
end
end
