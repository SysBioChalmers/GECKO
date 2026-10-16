function x = parseKcatNumber(texts)
% parseKcatNumber  Read numbers written in kcat source files.
%
% Accepts plain decimal numbers with an optional exponent (12, 0.5, 1e-05,
% 2.5E+3) and Inf; anything else, such as '5,5' from a spreadsheet with a
% decimal comma, becomes NaN instead of being misread. geckopy's
% read_kcat_source accepts the same text.
x = nan(numel(texts), 1);
ok = ~cellfun(@isempty, regexp(texts(:), '^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$|^[-+]?Inf$', 'once'));
x(ok) = str2double(texts(ok));
end
