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
%   model       an ecModel in GECKO 3 format (with ecModel.ec structure)
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
%   model       ecModel with updated model.ec.kcat and model.ec.source
%   rxnIdx      list of reaction indices (matching model.ec.rxns), to
%               indicate which kcat values have been changed.
% Usage:
%   [model, rxnIdx] = selectKcatValue(model,kcatList,criteria,overwrite)

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

% Remove zero kcat values. Only adjusting fields that are used later.
removeZero                      = kcatList.kcats == 0;
kcatList.kcats(removeZero)      = [];
kcatList.rxns(removeZero)       = [];

% Map to model.ec.rxns
[sanityCheck,idxInModel] = ismember(kcatList.rxns,model.ec.rxns);
if ~all(sanityCheck)
    error('Not all reactions in kcatList are found in model.ec.rxns')
end
% Make vector with single kcat value per reaction
idxInModelUnique = unique(idxInModel);
selectedKcats    = zeros(numel(idxInModelUnique),1);
selectedSource   = cell(numel(selectedKcats),1);
if ~isfield(kcatList,'kcatSource')
    kcatList.kcatSource = cell(numel(kcatList.kcats),1);
    kcatList.kcatSource(:) = {kcatList.source};
end
for i=1:numel(idxInModelUnique)
    ind = idxInModelUnique(i);
    idxMatch = find(idxInModel == ind);
    % Choose the maximum number
    switch criteria
        case 'max'
            [selectedKcats(i),j] = max(kcatList.kcats(idxMatch));
        case 'min'
            [selectedKcats(i),j] = min(kcatList.kcats(idxMatch));
        case 'median'
            [selectedKcats(i),j] = median(kcatList.kcats(idxMatch));
        case 'mean'
            [selectedKcats(i),j] = mean(kcatList.kcats(idxMatch));
        otherwise
            error('Invalid criteria specified')
    end
    selectedSource(i)    = kcatList.kcatSource(idxMatch(j));
end

% Populate model.ec.kcat
switch overwrite
    case 'true'
        model.ec.kcat(idxInModelUnique) = selectedKcats;
        model.ec.source(idxInModelUnique) = selectedSource;
    case 'false'
        emptyKcats = find(model.ec.kcat == 0);
        [idxInModelUnique,whickKcats] = intersect(idxInModelUnique,emptyKcats,'stable');
        model.ec.kcat(idxInModelUnique) = selectedKcats(whickKcats);
        model.ec.source(idxInModelUnique) = selectedSource(whickKcats);
        
    case 'ifHigher'
        higherKcats = model.ec.kcat(idxInModelUnique) < selectedKcats;
        selectedKcats(~higherKcats) = [];
        selectedSource(~higherKcats) = [];
        idxInModelUnique(~higherKcats) = [];
        model.ec.kcat(idxInModelUnique) = selectedKcats;
        model.ec.source(idxInModelUnique) = selectedSource;
    otherwise
        error('Invalid overwrite flag specified')
end
rxnIdx = idxInModelUnique;
end
