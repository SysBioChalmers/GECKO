function protData = loadProtData(replPerCond, varargin)
% loadProtData  Load absolute proteomics data and average over replicates.
%
% Loads absolute proteomics data (in mg/gDCW) and returns mean values
% across replicates for each condition in the data file. By default it also
% filters the data by various criteria, to remove uncertain data (see input
% parameters).
%
% Parameters
% ----------
% replPerCond : double
%     vector with number of replicates for each condition in the dataset.
%     Example: [3, 2] if first condition has triplicates and second
%     condition has duplicates.
%
% Name-Value Arguments
% --------------------
% protDataFile : char
%     path to file with proteomics data, where protein levels are in
%     mg/gDCW (default reads data/proteomics.tsv as specified in
%     modelAdapter). Alternatively, protDataFile can be a protData structure
%     that was previously made by loadProtData.
% filterData : logical
%     whether abundances should be filtered. If false, minVal, maxRSD,
%     maxMissing, cutLowest and addStdevs are not considered (default true).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% minVal : double
%     minimum mean protein measurement per condition; proteins with a lower
%     mean abundance are discarded (set to NaN) for that condition (default
%     0, i.e. no lower threshold).
% maxRSD : double
%     maximum relative standard deviation (RSD, the standard deviation
%     across replicates divided by their mean) allowed per condition;
%     proteins exceeding it are discarded (set to NaN) for that condition
%     (default 1, i.e. RSD up to 100%).
% maxMissing : double
%     minimum fraction of replicates within a condition that must have a
%     measured value; proteins with fewer present replicates than this
%     fraction requires are discarded (set to NaN) for that condition
%     (default 2/3, i.e. at least two of three replicates must have a
%     value). If conditions have different number of replicates (as
%     indicated in replPerCond), maxMissing can also be a vector of the
%     same length as replPerCond, with individual values per condition.
% cutLowest : double
%     percentage of lowest mean values per condition to be discarded (not
%     considering NaN values) (default 5).
% addStdevs : double
%     how many standard deviations should be added to the mean value of
%     each protein measurement across replicates, broadening the confidence
%     interval (default 1).
%
% Returns
% -------
% protData : struct
%     structure with (filtered) proteome data.
%
% Notes
% -----
% The protData structure has the following fields:
%
% - uniprotIDs : cell array with Uniprot IDs matching protData.abundances.
% - abundances : matrix of proteomics data, where each column contains mean abundances per condition.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     protData = loadProtData(replPerCond, protDataFile, filterData, modelAdapter, minVal, maxRSD, maxMissing, cutLowest, addStdevs);
%     protData = loadProtData(replPerCond, 'filterData', false);
%
% See also
% --------
% loadFluxData

p = parseGECKOargs(varargin, { ...
    'protDataFile', []; ...
    'filterData',   true; ...
    'modelAdapter', []; ...
    'minVal',       0; ...
    'maxRSD',       1; ...
    'maxMissing',   2/3; ...
    'cutLowest',    5; ...
    'addStdevs',    1});
protDataFile = p.protDataFile;
filterData   = p.filterData;
modelAdapter = p.modelAdapter;
minVal       = p.minVal;
maxRSD       = p.maxRSD;
maxMissing   = p.maxMissing;
cutLowest    = p.cutLowest;
addStdevs    = p.addStdevs;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(protDataFile)
    protDataFile = fullfile(params.path,'data','proteomics.tsv');
end

format = '%s';
for i=1:sum(replPerCond)
    format = [format ' %f'];
end
if ~isstruct(protDataFile)
    checkFileExistence(protDataFile);
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

