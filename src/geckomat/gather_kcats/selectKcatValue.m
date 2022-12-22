function [model, rxnIdx] = selectKcatValue(model,kcatList,criteria,overwrite)
% selectKcatValue
%   From a kcatList with predicted or suggested kcat values, where each
%   reaction may have multiple entries, one kcat value is selected and
%   written to model.ec.kcat. Zero values are discarded from the start. By
%   default, the maximum value is chosen, but alternatives are available.
%   The kcatList structure is an output of e.g. readDLKcatOutput,
%   readGotEnzymesOutput, readManualKcatList.
%
% Input:
%   model       an ec-model in RAVEN format
%   kcatList    structure array with separate entries for each kcat value
%               source      e.g. 'DLKcat' or 'gotenzymes'           
%               rxns        reaction identifiers, matching model.rxns
%               genes       gene identifiers, matching model.genes
%               substrate   substrates, matching model.mets
%               kcat        predicted kcat value in /sec
%   criteria    which kcat value should be selected if multiple values are
%               provided. Options: 'max', 'min', 'median', 'mean'. (Opt,
%               default 'max')
%   overwrite   whether existing kcat values should be overwritten.
%               Options: 'true', 'false', 'ifHigher'. The last option will
%               overwrite only if the new kcat value is higher. (Opt,
%               default 'true')
%
% Output:
%   model       ec-model with updated model.ec.kcat and model.ec.source
%   rxnIdx      list of reaction indices (matching model.ec.rxns), to
%               indicate which kcat values have been changed.
%   

if nargin < 4
    overwrite = 'true';
elseif islogical(overwrite)
    if overwrite
        overwrite = 'true';
    else
        overwrite = 'false';
    end
end
if nargin < 3
    criteria = 'max';
end

% Remove zero kcat values
removeZero                      = kcatList.kcats == 0;
kcatList.kcats(removeZero)      = [];
kcatList.rxns(removeZero)       = [];
kcatList.genes(removeZero)      = [];
kcatList.substrates(removeZero) = [];

% Map to model.ec.rxns
[sanityCheck,idxInModel] = ismember(kcatList.rxns,model.ec.rxns);
if ~all(sanityCheck)
    error('Not all reactions in kcatList are found in model.ec.rxns')
end
% Make vector with single kcat value per reaction
idxInModelUnique = unique(idxInModel);
selectedKcats    = zeros(numel(idxInModelUnique),1);
for i=1:numel(idxInModelUnique)
    ind = idxInModelUnique(i);
    % Choose the maximum number
    switch criteria
        case 'max'
            selectedKcats(i) = max(kcatList.kcats(idxInModel == ind));
        case 'min'
            selectedKcats(i) = min(kcatList.kcats(idxInModel == ind));
        case 'median'
            selectedKcats(i) = median(kcatList.kcats(idxInModel == ind));
        case 'mean'
            selectedKcats(i) = mean(kcatList.kcats(idxInModel == ind));
        otherwise
            error('Invalid criteria specified')
    end
end

% Populate model.ec.kcat
switch overwrite
    case 'true'
        model.ec.kcat(idxInModelUnique) = selectedKcats;
    case 'false'
        emptyKcats = find(model.ec.kcat == 0);
        [idxInModelUnique,whickKcats] = intersect(idxInModelUnique,emptyKcats,'stable');
        selectedKcats=selectedKcats(whickKcats);
        model.ec.kcat(idxInModelUnique) = selectedKcats;
    case 'ifHigher'
        higherKcats = model.ec.kcat(idxInModelUnique) < selectedKcats;
        selectedKcats(~higherKcats) = [];
        idxInModelUnique(~higherKcats) = [];
        model.ec.kcat(idxInModelUnique) = selectedKcats;
    otherwise
        error('Invalid overwrite flag specified')
end
model.ec.source(idxInModelUnique) = {kcatList.source};
rxnIdx = idxInModelUnique;
end
