function [baseRxns, enzymeSets] = getEcRxnKeys(model)
% getEcRxnKeys  Build-independent keys of the model.ec reactions.
%
% A kcat source file (see saveKcatSource) identifies a reaction by its
% conventional reaction and its enzyme set, not by model.ec.rxns, whose
% identifiers depend on how the ecModel was built. This function gives
% both keys for every model.ec row.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Returns
% -------
% baseRxns : cell
%     per model.ec.rxns entry, the conventional reaction: without a
%     trailing _EXP_<n> (full model) or a leading <nnn>_ (light model). A
%     reverse direction keeps its _REV.
% enzymeSets : cell
%     per model.ec.rxns entry, the UniProt accessions of its enzymes
%     (model.ec.rxnEnzMat), sorted and joined by ' + '; empty without
%     enzymes.
%
% Examples
% --------
%     [baseRxns, enzymeSets] = getEcRxnKeys(ecModel);
%
% See also
% --------
% saveKcatSource, loadKcatSource

if model.ec.geckoLight
    baseRxns = regexprep(model.ec.rxns, '^\d+_', '');
else
    baseRxns = regexprep(model.ec.rxns, '_EXP_\d+$', '');
end
baseRxns = baseRxns(:);

n = numel(model.ec.rxns);
enzymeSets = repmat({''}, n, 1);
[rows, cols] = find(model.ec.rxnEnzMat);
if ~isempty(rows)
    % Enzymes per reaction in one pass, rather than a search per reaction.
    perRxn = accumarray(rows(:), cols(:), [n 1], @(c) {c});
    for i = find(~cellfun(@isempty, perRxn))'
        enzymes = unique(model.ec.enzymes(perRxn{i}));
        enzymeSets{i} = strjoin(enzymes(:)', ' + ');
    end
end
end
