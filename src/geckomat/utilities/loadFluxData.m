function fluxData = loadFluxData(fluxDataFile, modelAdapter)
% loadFluxData
%   Function that loads total protein measurements and flux data (exchange
%   fluxes for carbon source,  O2, CO2, etc.)
%
% Input:
%   fluxDataFile    path to file with flux data. (Optional, default
%                   reads data/fluxData.tsv as specified in modelAdapter)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use
%                   the default model adapter)
%
% Output:
%   fluxData                structure with flux data
%       conds               sampling condition
%       Ptot                total protein (g/gDCW)
%       grRate              growth rate (1/h)
%       exchFluxes          exchange fluxes (mmol/gDCWh)
%       exchMets            exchanged metabolites, matching exchFluxes
%       exchRxnIDs          exchange reaction IDs, matching exchMets
%       bayesianRMSEweight  if column existed in fluxData.tsv, weights for
%                           RMSE calculation in Bayesian kcat tuning
%       source              if column existed in fluxData.tsv, description
%                           where the data comes from
%
% Usage:
%   fluxData = loadFluxData(fluxDataFile, modelAdapter)

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 1 || isempty(fluxDataFile)
    fluxDataFile = fullfile(params.path,'data','fluxData.tsv');
end

%Load total protein content and flux data
fID       = fopen(fluxDataFile);
formatStr = '%s';
data      = textscan(fID,formatStr,'Delimiter','\n');
fclose(fID);
fluxDataRaw  = [];
for i=1:length(data{1})
    row      = data{1}(i);
    row      = strsplit(row{1},'\t');
    row      = row(1:end);
    fluxDataRaw = [fluxDataRaw; row]; 
end
if size(fluxDataRaw,2) == 0
    error('The input file %s is likely not tab-delimited.', fluxDataFile)
end
fluxData            = [];

%Find additional fields
extraFieldNames     = {'bayesianRMSEweight','source'};
extraFieldNumeric   = [true,false];      
[~, extraFields]    = ismember(extraFieldNames,fluxDataRaw(1,:));
for i=1:numel(extraFields)
    if extraFields(1) ~= 0
        if extraFieldNumeric(i)
            fluxData.(extraFieldNames{i}) = str2double(fluxDataRaw(2:end,extraFields(i)));
        else
            fluxData.(extraFieldNames{i}) = fluxDataRaw(2:end,extraFields(i));
        end
    end
end
fluxDataRaw(:,extraFields) = [];

%Extract observed byProduct names from file
exchRxns = fluxDataRaw(1,4:end);
exchMets = strtrim(regexprep(exchRxns,'(.*)\(.*\)$','$1'));
exchRxns = regexprep(exchRxns,'.*\((.*)\)$','$1');

fluxData.conds      = fluxDataRaw(2:end,1);
fluxData.Ptot       = str2double(fluxDataRaw(2:end,2));
fluxData.grRate     = str2double(fluxDataRaw(2:end,3));
fluxData.exchFluxes = str2double(fluxDataRaw(2:end,4:end));
%fluxData.exchFluxes(isnan(fluxData.exchFluxes)) = 0;
fluxData.exchMets   = exchMets;
fluxData.exchRxnIDs = exchRxns;
end
