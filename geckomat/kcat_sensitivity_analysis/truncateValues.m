%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% table = truncateValues(table,nCols)
%
%
% Benjamin Sanchez    Last edited: 2018-05-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table = truncateValues(table,nCols)

[m,n] = size(table);
for i = 1:m
    for j = (n-nCols+1):n
        orderMagn  = max([ceil(log10(abs(table{i,j}))),0]);
        table{i,j} = round(table{i,j},6-orderMagn);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
