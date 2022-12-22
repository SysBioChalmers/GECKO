function [pIDs,abundances] = filter_ProtData(uniprots,data,flexFactor,medianF,minVal)
% filter_ProtData
%
% Function that takes a proteomics dataset for a given condition and
% filters out those proteins with empty IDs, absent in more than 1/3 of the
% biological replicates and with a RSD>1
%
%   uniprots    (cell) Uniprot IDs for a proteomics dataset
%   data        (matrix) Numerical values for protein abundances (absolute
%               or relative [number of proteins x number of samples]
%   flexFactor  (num) Fraction of the standard deviation that should be
%               added to the mean value of each protein measurement across 
%               replicates (confidence interval).
%   medianF     True if the RSD should be calculated with respect to the
%               median value of a protein abundance (stringent filter) or 
%               False if the mean value is used instead.
%   minVal      Detection limit for proteins (default = 0).
%
% Created.       Ivan Domenzain 2019-02-13
% Last modified. Ivan Domenzain 2019-04-12

if nargin<5
    minVal = 0;
    if nargin<4
        medianF = false;
        if nargin<3
            flexFactor = 1;
        end
    end
end
pIDs       = [];
abundances = [];
[m,n]      = size(data);
for i=1:m
    %Check for non empty ID
    if ~isempty(uniprots{i})
        row = data(i,:);
        %Proteins should be present in, at least, 2/3 of the replicates
        if sum(row>0)>=(2/3)*n
            %variability filter
            if medianF
                metric = std(row)/median(row);
            else
                metric = std(row)/mean(row);
            end
            value = mean(row)+(flexFactor*std(row));
            if metric<1 && value>minVal
                abundances = [abundances; value];
                pIDs       = [pIDs;uniprots(i)];
            end
        end
    end
end
end