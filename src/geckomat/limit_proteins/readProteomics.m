function [model, protData] = readProteomics(model, condition, protData, modelAdapter)
% readProteomics
%   Reads absolute proteomics data, from either a file or protData
%   structure. The proteins should be annotated with Uniprot IDs, while the
%   protein levels should be given as mmol/gDCW. Any existing protein
%   concentrations in model.ec.concs are discarded. If no data is provided
%   a particular protein, its level is NaN.
%
% Input:
%   model           an ec-model
%   condition       a condition to be set as protein levels. It must be a
%                   number (Optional, default = 1)
%   protData        either string where to find the input file, or a
%                   structure with the following fields:
%                   protData.uniprot  Uniprot identifiers
%                   protData.level    Protein levels in mmol/gDCW. 
%                   If it is a file multiple conditions can be load.
%                   However, only the condition defined will be used in
%                   model.ec.concs
%   modelAdapter    a loaded model adapter (Optional, will otherwise use
%                   the default model adapter)
% Output:
%   model           an ec-model where model.ec.concs is populated with
%                   protein concentrations.
%
% Note: to also constrain the model with the content of model.ec.concs, you
% should run constrainProtConcs.

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 3 || isempty(protData)
    protData = fullfile(params.path,'data','proteomics.tsv');
end

if nargin < 2
    condition = 1;
end

% Chack if condition is a numeric type and define the set data to used
if isnumeric(condition)
    n = 1 + condition;
    % Repeat the predefined format as many conditions
    format = join(repmat({'%s'},1,n),' ');
    % If there is more than n columns, skip them
    format = [format{:} '%*[^\n]']; 
else
    error('Condition must be a number')
end

switch class(protData)
    case {'char','string'}
        fID = fopen(fullfile(protData),'r');
        fileContent = textscan(fID, format, 'HeaderLines', 1);
        fclose(fID);
        uniprot = fileContent{1};
        level   = str2double(fileContent{n});
    case 'struct'
        uniprot = protData.uniprot;
        level   = protData.level;
    otherwise
        error(['protData should either be a structure or string, not a ' class(protData)])
end

%Redefine an empty model.ec.concs vector
model.ec.concs=nan(numel(model.ec.enzymes),1);

[a,b] = ismember(uniprot, model.ec.enzymes);
model.ec.concs(b(a)) = level(a);

clear protData
protData.uniprot = uniprot;
protData.level = level;
end
