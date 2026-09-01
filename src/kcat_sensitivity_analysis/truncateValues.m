function table = truncateValues(table,cols)
% truncateValues  Round values in selected columns to six significant digits.
%
% Rounds each value in the given columns of a cell array or table to six
% significant digits, based on the value's order of magnitude.
%
% Parameters
% ----------
% table : cell or table
%     a cell array or table where some columns have values that should be
%     truncated.
% cols : double
%     index or indices of columns with values to be truncated.
%
% Returns
% -------
% table : cell or table
%     the input with the selected columns rounded.

[m,n] = size(table);
for i = 1:m
    for j = cols
        orderMagn  = max([ceil(log10(abs(table{i,j}))),0]);
        table{i,j} = round(table{i,j},6-orderMagn);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
