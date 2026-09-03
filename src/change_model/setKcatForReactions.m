function ecModel = setKcatForReactions(ecModel,rxnIds,kcat)
% setKcatForReactions  Change the kcat value for selected reactions.
%
% Changes the kcat value in ecModel.ec.kcat for selected reactions.
% applyKcatConstraints needs to be run afterwards to transfer the kcat
% values into the S-matrix.
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% rxnIds : char or cell
%     reaction identifier matching ecModel.ec.rxns. If the _EXP_ suffix is
%     not included, and there are multiple expanded (isozymic) reactions,
%     then all of those will have their kcat changed. If rxnIds includes a
%     _EXP_ suffix, then only that specific reaction will have its kcat
%     changed. If multiple rxnIds are provided as a cell array, then the
%     above applies to each rxnIds individually.
% kcat : double
%     the new kcat value.
%
% Returns
% -------
% ecModel : struct
%     ecModel where selected kcat values in ecModel.ec.kcat are changed, but
%     not yet applied to the S-matrix (will require to run
%     applyKcatConstraints). ecModel.ec.source for the changed reactions
%     will read 'setKcatForReactions'.
%
% Examples
% --------
%     ecModel = setKcatForReactions(ecModel, rxnIds, kcat);
%
% See also
% --------
% applyKcatConstraints, getKcatAcrossIsozymes
rxnIds = convertCharArray(rxnIds);

hasExp       = ~cellfun(@isempty,regexp(rxnIds,'_EXP_\d+$'));
nonExpRxns   = regexprep(ecModel.ec.rxns,'_EXP_\d+$','');
rxnsToChange = [];
for i=1:numel(hasExp)
    if hasExp(i) == 1
        rxnsToChange = [rxnsToChange; find(strcmpi(ecModel.ec.rxns,rxnIds{i}))];
    else
        nonExpRxn    = regexprep(rxnIds(i),'_EXP_\d+$','');
        rxnsToChange = [rxnsToChange; find(strcmpi(nonExpRxns,nonExpRxn))];
    end
end
if isscalar(rxnsToChange)
    if length(kcat) ~= 1
        error('Found one reaction whose kcat should change, you should provide one kcat value only.')
    end
else
    if isscalar(kcat)
        % Is fine, all reactions get the same kcat
    elseif length(kcat) ~= length(rxnsToChange)
        error('Found %d reactions whose kcat should change, the new kcat should be either a single value, or a vector of length %d.', length(rxnsToChange), length(rxnsToChange))
    end
end
ecModel.ec.kcat(rxnsToChange)   = kcat;
ecModel.ec.source(rxnsToChange) = {'setKcatForReactions'};
end
