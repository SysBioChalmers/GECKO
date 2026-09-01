function ecModel = applyECcodes(ecModel, varargin)
% applyECcodes  Copy EC numbers from ecModel.ec.eccodes to ecModel.eccodes.
%
% Copies the EC numbers stored in the ecModel.ec structure across to the
% top-level ecModel.eccodes field.
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% overwrite : logical
%     whether existing (non-empty) ecModel.eccodes entries should be
%     overwritten by a non-empty ecModel.ec.eccodes entry. Empty
%     ecModel.eccodes entries are always overwritten if possible (default
%     false). An empty ecModel.ec.eccodes source entry is never copied, so
%     it can never erase an existing annotation.
%
% Returns
% -------
% ecModel : struct
%     ecModel with populated model.eccodes.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     ecModel = applyECcodes(ecModel, true);
%     ecModel = applyECcodes(ecModel, 'overwrite', true);

p = parseGECKOargs(varargin, { ...
    'overwrite', false});
overwrite = p.overwrite;

if ~isfield(ecModel.ec,'eccodes')
    error('ecModel.ec.eccodes does not exist')
end

[a,b] = ismember(ecModel.rxns,ecModel.ec.rxns);
if ~isfield(ecModel,'eccodes')
    ecModel.eccodes = cell(numel(ecModel.rxns),1);
end

if overwrite
    % Never overwrite an existing eccodes entry with an empty source: an
    % unset ec.eccodes value means "no info to copy", not "clear the
    % annotation".
    nonEmptySource = false(size(a));
    nonEmptySource(a) = ~cellfun(@isempty, ecModel.ec.eccodes(b(a)));
    a = a & nonEmptySource;
    ecModel.eccodes(a) = ecModel.ec.eccodes(b(a));
else % Only replace emptyEcCodes
    emptyEcCodes = cellfun(@isempty, ecModel.eccodes);
    a = a & emptyEcCodes;
    ecModel.eccodes(a) = ecModel.ec.eccodes(b(a));
end
end
