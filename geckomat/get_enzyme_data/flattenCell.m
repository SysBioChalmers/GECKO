function Cflat = flattenCell(C,strFlag)
%FLATTENCELL  Flatten a nested column cell array into a matrix cell array.
%
% CFLAT = FLATTENCELL(C) takes a column cell array in which one or more
% entries is a nested cell array, and flattens it into a 2D matrix cell
% array, where the nested entries are spread across new columns.
%
% CFLAT = FLATTENCELL(C,STRFLAG) if STRFLAG is TRUE, empty entries in the
% resulting CFLAT will be replaced with empty strings {''}. Default = FALSE
%
% NOTE: This function is scheduled to be added to a new release of the 
%       RAVEN toolbox, i.e. can be removed from GECKO after that:
%       https://github.com/SysBioChalmers/RAVEN/blob/d992169e9da44ed9c441753bcf77cc83f521e953/external/MetaNetX/flattenCell.m
%
% Jonathan Robinson, 2018-03-07

if nargin < 2
    strFlag = false;
end

% determine which entries are cells
cells = cellfun(@iscell,C);

% determine number of elements in each nested cell
cellsizes = cellfun(@numel,C);
cellsizes(~cells) = 1;  % ignore non-cell entries

% flatten single-entry cells
Cflat = C;
Cflat(cells & (cellsizes == 1)) = cellfun(@(x) x{1},Cflat(cells & (cellsizes == 1)),'UniformOutput',false);

% iterate through multi-entry cells
multiCells = find(cellsizes > 1);
for i = 1:length(multiCells)
    cellContents = Cflat{multiCells(i)};
    Cflat(multiCells(i),1:length(cellContents)) = cellContents;
end

% change empty elements to strings, if specified
if ( strFlag )
    Cflat(cellfun(@isempty,Cflat)) = {''};
end
