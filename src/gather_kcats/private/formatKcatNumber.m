function t = formatKcatNumber(x)
% formatKcatNumber  Number text of kcat source files.
%
% Whole numbers as integers, others in the shortest decimal form that reads
% back as the same double; 'Inf' or '-Inf' for infinity, empty for NaN.
% geckopy's format_number writes the same text.
if isnan(x)
    t = '';
elseif isinf(x)
    t = 'Inf';
    if x < 0
        t = '-Inf';
    end
elseif x == round(x) && abs(x) < 1e15
    t = sprintf('%d', x);
else
    for digits = 1:17
        t = sprintf('%.*g', digits, x);
        if str2double(t) == x
            break
        end
    end
end
end
