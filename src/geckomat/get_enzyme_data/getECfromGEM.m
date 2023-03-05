function model = getECfromGEM(model, ecRxns)
% getECfromGEM
%   Use the model.eccodes to populates the model.ec.eccodes field.
%
% Input:
%   model           ec-model in GECKO 3 format
%   ecRxns          logical of length model.ec.rxns that specifies for
%                   which reactions the existing model.ec.eccodes entry
%                   should be kept and not modified by this function
%                   (optional, by default all model.ec.eccodes entries
%                   are populated by this function)
%
% Output:
%   model           ec-model with populated model.ec.eccodes

if ~isfield(model,'eccodes')
    error('The model has no model.eccodes field.')
end

%Need to remove the prefix of GECKO light rxn names in the ec structure
if ~model.ec.geckoLight
    rxnNames = model.ec.rxns;
else
    rxnNames = extractAfter(model.ec.rxns, 4);
end

rxnIdxs = getIndexes(model,rxnNames,'rxns');

% Check if eccodes are valid
eccodes = model.eccodes;
noEcCodes = cellfun(@isempty, eccodes);
eccodes = eccodes(~noEcCodes);
rxns = model.rxns(~noEcCodes);
invalidEc = regexprep(eccodes,'(\d\.(\w|-)+\.(\w|-)+\.(\w|-)+)(;\w+\.(\w|-)+\.(\w|-)+\.(\w|-)+)*(.*)','$3');
invalidEc = ~cellfun(@isempty,invalidEc);
if any(invalidEc)
    fprintf(['The following reactions have invalid/incomplete EC numbers, correct them in model.eccodes and rerun getECfromGEM:\n%s\n'], ...
         strjoin(rxns(invalidEc),'\n'))
    return
end
if nargin<2 || all(ecRxns)
    model.ec.eccodes = model.eccodes(rxnIdxs);
else
    if ~isfield(model.ec,'eccodes')
        model.ec.eccodes(1:numel(model.ec.rxns),1) = {''};
    end
    model.ec.eccodes(ecRxns) = model.eccodes(rxnIdxs(ecRxns));
end
end
