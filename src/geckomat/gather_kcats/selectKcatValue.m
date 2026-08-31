function [model, rxnIdx] = selectKcatValue(model,kcatList,varargin)
% selectKcatValue  Select one kcat value per reaction and write to model.ec.
%
% From a kcatList with predicted or suggested kcat values, where each
% reaction may have multiple entries, one kcat value is selected and written
% to model.ec.kcat. Zero values are discarded from the start. By default,
% the maximum value is chosen, but alternatives are available. The kcatList
% structure is an output of e.g. readDLKcatOutput, readGotEnzymesOutput,
% readManualKcatList.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% kcatList : struct
%     structure array with separate entries for each kcat value, with the
%     following fields:
%
%     - source : e.g. 'DLKcat' or 'gotenzymes'.
%     - rxns : reaction identifiers, matching model.ec.rxns.
%     - genes : gene identifiers, matching model.ec.genes.
%     - substrates : substrate names, matching model.metNames.
%     - kcats : predicted kcat value in /sec.
%
% Name-Value Arguments
% --------------------
% criteria : char
%     which kcat value should be selected if multiple values are provided.
%     Options: 'max', 'min', 'median', 'mean' (default 'max'). The written
%     model.ec.source is taken from the matching entry for 'max'/'min'; for
%     'median'/'mean', where the selected value need not equal any single
%     input, it is taken from the first entry instead (arbitrary but
%     deterministic, since an aggregate has no single "winning" source).
% overwrite : char
%     whether existing kcat values should be overwritten. Options: 'true',
%     'false', 'ifHigher'. The last option will overwrite only if the new
%     kcat value is higher (default 'true').
%
% Returns
% -------
% model : struct
%     ecModel with updated model.ec.kcat and model.ec.source.
% rxnIdx : double
%     list of reaction indices (matching model.ec.rxns), to indicate which
%     kcat values have been changed.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, rxnIdx] = selectKcatValue(model, kcatList);
%     [model, rxnIdx] = selectKcatValue(model, kcatList, 'criteria', 'mean');
%
% See also
% --------
% readDLKcatOutput, mergeDLKcatAndFuzzyKcats, applyKcatConstraints

p = parseGECKOargs(varargin, { ...
    'criteria',  []; ...
    'overwrite', []});
criteria  = p.criteria;
overwrite = p.overwrite;

if isempty(overwrite)
    overwrite = 'true';
elseif islogical(overwrite)
    if overwrite
        overwrite = 'true';
    else
        overwrite = 'false';
    end
end
if isempty(criteria)
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
            selectedKcats(i) = median(kcatList.kcats(idxMatch));
            j = 1;
        case 'mean'
            selectedKcats(i) = mean(kcatList.kcats(idxMatch));
            j = 1;
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
