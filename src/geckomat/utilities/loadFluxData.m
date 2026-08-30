function fluxData = loadFluxData(varargin)
% loadFluxData  Load total protein measurements and exchange flux data.
%
% Loads total protein content and exchange flux data (e.g. carbon source,
% O2, CO2) from a tab-delimited file.
%
% Name-Value Arguments
% --------------------
% fluxDataFile : char
%     path to file with flux data (default reads data/fluxData.tsv as
%     specified in modelAdapter).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% fluxData : struct
%     structure with flux data.
%
% Notes
% -----
% The fluxData structure has the following fields:
%
% - conds : sampling condition.
% - Ptot : total protein (g/gDCW).
% - grRate : growth rate (1/h).
% - exchFluxes : exchange fluxes (mmol/gDCWh).
% - exchMets : exchanged metabolites, matching exchFluxes.
% - exchRxnIDs : exchange reaction IDs, matching exchMets.
% - bayesianRMSEweight : if column existed in fluxData.tsv, weights for RMSE calculation in Bayesian kcat tuning.
% - source : if column existed in fluxData.tsv, description where the data comes from.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     fluxData = loadFluxData(fluxDataFile, modelAdapter);
%     fluxData = loadFluxData('fluxDataFile', fluxDataFile);
%
% See also
% --------
% loadProtData

p = parseGECKOargs(varargin, { ...
    'fluxDataFile', []; ...
    'modelAdapter', []});
fluxDataFile = p.fluxDataFile;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(fluxDataFile)
    fluxDataFile = fullfile(params.path,'data','fluxData.tsv');
end

%Load total protein content and flux data
try
    fID       = fopen(fluxDataFile);
    formatStr = '%s';
    data      = textscan(fID,formatStr,'Delimiter','\n');
    fclose(fID);
catch
    warning('Failed to load %s.',fluxDataFile);
    fluxData = [];
    return
end
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
[logicalExtra, extraFields]    = ismember(extraFieldNames,fluxDataRaw(1,:));
extraFields(~logicalExtra) = [];
extraFieldNames(~logicalExtra) = [];
extraFieldNumeric(~logicalExtra) = [];
if any(extraFields)
    for i=1:numel(extraFields)
        if extraFieldNumeric(i)
            fluxData.(extraFieldNames{i}) = str2double(fluxDataRaw(2:end,extraFields(i)));
        else
            fluxData.(extraFieldNames{i}) = fluxDataRaw(2:end,extraFields(i));
        end
    end
    fluxDataRaw(:,extraFields) = [];
end

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
