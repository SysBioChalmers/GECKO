function protData = loadProtData(replPerCond, protDataFile, filterData, modelAdapter, minVal, maxRSD, maxMissing, cutLowest, addStdevs)
% loadProtData
%   Function that loads absolute proteomics data (in mg/gDCW) and returns
%   mean values across replicates for each condition in the data file. By
%   default it also filters the data by various criteria, to remove
%   uncertain data (see input parameters).
%
% Input:
%   replPerCond     vector with number of replicates for each condition in
%                   the dataset. Example: [3, 2] if first conditions has
%                   triplicates and second condition has duplicates.
%   protDataFile    path to file with proteomics data, where protein levels
%                   are in mg/gDCW (Optional, default reads 
%                   data/proteomics.tsv as specified in modelAdapter)
%                   Alternatively, protDataFile can be a protData structure
%                   that was previously made by loadProtdata.
%   filterData      logical whether abundances should be filtered. If
%                   false, minVal, maxRSD, maxMissing and addStdevs are not
%                   considered. (Optional, default true)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use
%                   the default model adapter)
%   minVal          threshold of mean protein measurement per condition.
%                   (Optional, default = 0)
%   maxRSD          maximum relative standard per condition. (Optional,
%                   default = 1)
%   maxMissing      ratio of replicates for which a protein level might be
%                   missing. (Optional, default = 1/3 (or 1/2 if number of
%                   replicates = 2))
%                   If conditions have different number of replicates (as
%                   indicated in replPerCond), maxMissing can also be a
%                   vector of the same length as replPerCond, with
%                   individual maxMissing parameters for each replicate.
%   cutLowest       percentage of lowest mean values per condition to be
%                   discared (not considering NaN values). (Optional, default 5)
%   addStdevs       how many standard deviations should be added to the mean
%                   value of each protein measurement across replicates,
%                   broadening the confidence interval. (Optional,
%                   default = 1)
%
% Output:
%   protData        structure with (filtered) proteome data
%                   uniprotIDs  cell arrray with Uniprot IDs matching
%                               protData.abundances
%                   abundances  matrix of proteomics data, where each
%                               column contains mean abundances per
%                               condition
%
% Usage:
%   protData = loadProtData(replPerCond, protDataFile, filterData, modelAdapter, minVal, maxRSD, maxMissing, cutLowest, addStdevs)

if nargin < 8 || isempty(addStdevs)
    addStdevs = 1;
end

if nargin < 7 || isempty(maxMissing)
    maxMissing = 2/3;
end

if nargin < 6 || isempty(maxRSD)
    maxRSD = 1;
end

if nargin < 5 || isempty(minVal)
    minVal = 0;
end

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 3 || isempty(filterData)
    filterData = true;
end

if nargin < 2 || isempty(protDataFile)
    protDataFile = fullfile(params.path,'data','proteomics.tsv');
end

format = '%s';
for i=1:sum(replPerCond)
    format = [format ' %f'];
end
if ~isstruct(protDataFile)
    fID         = fopen(protDataFile);
    protDataRaw = textscan(fID,format,'Delimiter','\t','HeaderLines',1,'TreatAsEmpty',{'NA','na','NaN','#VALUE!'});
    uniprotIDs  = protDataRaw{1};
    abundances  = cell2mat(protDataRaw(2:end));
    fclose(fID);
else
    uniprotIDs  = protDataFile.uniprotIDs;
    abundances  = protDataFile.abundances;
end

%Remove entriew without ID
remData = cellfun(@isempty,uniprotIDs);
uniprotIDs(remData,:) = [];
abundances(remData,:) = [];
m                     = size(abundances,1);
filtAbund             = nan(m,numel(replPerCond));

if filterData
    for i=1:numel(replPerCond)
        condAbund    = abundances(:,1:replPerCond(i));
        if i<numel(replPerCond)
            abundances   = abundances(:,replPerCond(i)+1:end);
        end
        % First filter maxMissing
        if size(condAbund,2) > 1
            if numel(maxMissing)>1
                maxMisRepl = maxMissing(i);
            else
                maxMisRepl = maxMissing;
            end
            remData = sum(condAbund>0,2)<maxMisRepl*size(condAbund,2);
            condAbund(remData,:) = nan;
        end
        % Filter by RSD
        remData = (std(condAbund,0,2,'omitnan')./mean(condAbund,2,'omitnan'))>maxRSD;
        condAbund(remData) = nan;
        % Add stdevs
        condAbund = mean(condAbund,2,'omitnan')+(addStdevs*std(condAbund,0,2,'omitnan'));
        % Filter by minVal
        remData = mean(condAbund,2,'omitnan')<minVal;
        condAbund(remData) = nan;
        % Remove bottom 5%
        numData  = find(~isnan(condAbund));
        [~,sortData] = sort(condAbund);
        lowCutoff = floor(numel(numData)*0.05);
        condAbund(sortData(1:lowCutoff)) = nan;
        % Combine conditions
        filtAbund(:,i) = condAbund;
    end
else
    for i=1:numel(replPerCond)
        condAbund    = abundances(:,1:replPerCond(i));
        if i<numel(replPerCond)
            abundances = abundances(:,replPerCond(i)+1:end);
        end
        filtAbund(:,i) = mean(condAbund,2,'omitnan');
    end
end
notAllNan = logical(sum(~isnan(filtAbund),2));
protData.abundances = filtAbund(notAllNan,:);
protData.uniprotIDs = uniprotIDs(notAllNan);
end

