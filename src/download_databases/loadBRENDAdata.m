%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KCATcell, SAcell] = loadBRENDAdata(varargin)
% loadBRENDAdata  Load kcat and specific activity data from BRENDA files.
%
% Reads the BRENDA data files (kcat.tsv, sa.tsv and mw.tsv) from the BRENDA
% database folder defined by the model adapter, and returns the kcat values
% and specific activities (the latter converted to kcat values using the
% molecular weights).
%
% kcat.tsv and sa.tsv each carry both a max and a median aggregate per
% (EC, substrate, organism) triple; kcatAggregation selects which of the
% two is read. Each file starts with a `#`-prefixed release-version line
% followed by a tab-delimited column header, both skipped on read.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% kcatAggregation : char
%     which aggregate to read from kcat.tsv and sa.tsv, 'max' or
%     'median' (default: params.kcatAggregation from the model adapter,
%     or 'max' if the adapter does not set it). A triple's max is the
%     highest turnover reported for it; its median is the middle one.
%
% Returns
% -------
% KCATcell : cell
%     kcat data extracted from BRENDA.
% SAcell : cell
%     specific activity data extracted from BRENDA.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [KCATcell, SAcell] = loadBRENDAdata();
%     [KCATcell, SAcell] = loadBRENDAdata('modelAdapter', adapter);
%     [KCATcell, SAcell] = loadBRENDAdata('kcatAggregation', 'median');
%
% See also
% --------
% loadDatabases, getECfromDatabase

p = parseGECKOargs(varargin, { ...
    'modelAdapter',    [], @(x) assert(isempty(x) || isa(x,'ModelAdapter'), ...
                                'loadBRENDAdata:invalidModelAdapter', ...
                                'modelAdapter must be a ModelAdapter.'); ...
    'kcatAggregation', [], @(x) assert(isempty(x) || any(strcmpi(x, ...
                                {'max','median'})), ...
                                'loadBRENDAdata:invalidKcatAggregation', ...
                                'kcatAggregation must be ''max'' or ''median''.')});
modelAdapter    = p.modelAdapter;
kcatAggregation = p.kcatAggregation;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

kcatAggregation = resolveKcatAggregation(kcatAggregation, modelAdapter);
%Column 4 holds the per-triple maximum and column 5 the median.
valueCol = 4 + strcmpi(kcatAggregation,'median');

basePath      = modelAdapter.getBrendaDBFolder();
KCAT_file      = fullfile(basePath,'kcat.tsv');
SA_file        = fullfile(basePath,'sa.tsv');
MW_file        = fullfile(basePath,'mw.tsv');

%Extract BRENDA DATA from files information. kcat.tsv/sa.tsv are seven
%columns (ec_code, substrate, organism, value_max, value_median, n,
%references); valueCol picks max or median. mw.tsv is six columns (no
%substrate-level aggregation choice): ec_code, substrate, organism, value,
%n, references.
KCATcell      = openDataFile(KCAT_file,1,'%q %q %q %f %f %f %q',valueCol);
scalingFactor = 1/60;    %[umol/min/mg] -> [mmol/s/g]
SA            = openDataFile(SA_file,scalingFactor,'%q %q %q %f %f %f %q',valueCol);
scalingFactor = 1/1000;  %[g/mol] -> [g/mmol]
MW            = openDataFile(MW_file,scalingFactor,'%q %q %q %f %f %q',4);

for i=1:4
    SAcell{i} = [];
end
previousEC = []; EC_indexes = [];

%build an index on MW{1} to speed things up a bit
%first just extract the genus (i.e. the first part of the name)
MWECNum = upper(unique(MW{1}));
MWECNumIndices = cell(length(MWECNum),1);
MWECNumHashMap = containers.Map(MWECNum,1:length(MWECNum));
for i = 1:length(MW{1})
    matchInd = cell2mat(values(MWECNumHashMap, MW{1}(i)));
    MWECNumIndices{matchInd} = [MWECNumIndices{matchInd};i];
end


for i=1:length(SA{1})
    %Gets the indexes of the EC repetitions in the MW cell for every
    %new (different) EC
    if ~strcmpi(SA{1}(i), previousEC)
        key = upper(SA{1}(i));
        if isKey(MWECNumHashMap,key) %annoyingly, this seems to be needed
            matchInd = cell2mat(values(MWECNumHashMap,key));
            EC_indexes = MWECNumIndices{matchInd};
        else
            EC_indexes = [];
        end
    end
    mwEC{1} = MW{3}(EC_indexes); mwEC{2} = MW{4}(EC_indexes);
    % just looks for the first match because just the maximal value for
    % each EC# / Orgaism is reported on the file
    org_index = find(strcmpi(SA{3}(i),mwEC{1}),1);
    if ~isempty(org_index)
        SAcell{1} = [SAcell{1};SA{1}(i)];
        SAcell{2} = [SAcell{2};SA{3}(i)];
        SAcell{3} = [SAcell{3}; SA{4}(i)*mwEC{2}(org_index)]; %[1/s]
        SAcell{4} = [SAcell{4}; mwEC{2}(org_index)];
    end
    previousEC = SA{1}(i);
end

function data_cell = openDataFile(fileName,scalingFactor,formatSpec,valueCol)
fID          = fopen(fileName);
raw          = textscan(fID,formatSpec,'delimiter','\t','HeaderLines',2);
fclose(fID);
data_cell    = {raw{1}, raw{2}, raw{3}, raw{valueCol}*scalingFactor};
end
end
