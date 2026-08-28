function [model, invalidEC, invalidECpos] = getECfromGEM(model, varargin)
% getECfromGEM  Populate model.ec.eccodes from the model.eccodes field.
%
% Uses the model.eccodes to populate the model.ec.eccodes field. EC numbers
% must be formatted as four numbers separated by periods, possibly with
% trailing wildcards; incorrectly formatted EC numbers are skipped.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% ecRxns : logical
%     of length model.ec.rxns that specifies for which reactions the
%     existing model.ec.eccodes entry should be kept and not modified by
%     this function (by default all model.ec.eccodes entries are populated
%     by this function).
%
% Returns
% -------
% model : struct
%     ecModel with populated model.ec.eccodes.
% invalidEC : cell
%     incorrectly formatted EC numbers.
% invalidECpos : double
%     position of invalidEC in model.eccodes.
%
% Notes
% -----
% Valid EC numbers are e.g. 1.2.3.4 or 1.2.3.-, while invalid EC numbers are
% 1.2.3 or 1_2_3_4. Multiple EC numbers are separated by ; for instance
% 1.2.3.4;1.2.3.5 not 1.2.3.4|1.2.3.5.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, invalidEC, invalidECpos] = getECfromGEM(model, ecRxns);
%     [model, invalidEC, invalidECpos] = getECfromGEM(model, 'ecRxns', ecRxns);
%
% See also
% --------
% getECfromDatabase, copyECtoGEM

p = parseGECKOargs(varargin, { ...
    'ecRxns', []});
ecRxns = p.ecRxns;

if ~isfield(model,'eccodes')
    error('The model has no model.eccodes field.')
end

%Need to remove the prefix of GECKO light rxn names in the ec structure
if ~model.ec.geckoLight
    rxnNames = model.ec.rxns;
else
    rxnNames = extractAfter(model.ec.rxns, 4);
end

[~,rxnIdxs] = ismember(rxnNames,model.rxns);

% Check if eccodes are valid. Each of the four dot-separated components must
% be either a run of digits or the wildcard '-'; \w (the previous class) also
% matches letters and underscores, so malformed codes such as '1.1.1.n12' or
% '1_2_3_4' were incorrectly accepted.
eccodes = model.eccodes;
invalidEC = regexprep(eccodes,'(\d+\.(\d+|-)\.(\d+|-)\.(\d+|-))(;\d+\.(\d+|-)\.(\d+|-)\.(\d+|-))*(.*)','$3');
invalidEC = ~cellfun(@isempty,invalidEC);
invalidECpos = find(invalidEC);
if any(invalidECpos)
    invalidEC = model.eccodes(invalidEC);
    if nargout<2
        fprintf('Skipped incorrectly formatted EC numbers, rerun getECfromGEM with all outputs to get a list.\n')
    else
        fprintf('Skipped incorrectly formatted EC numbers.\n')
    end
    eccodes(invalidECpos)={''};
else
    invalidEC = [];
end
if isempty(ecRxns) || all(ecRxns)
    model.ec.eccodes = eccodes(rxnIdxs);
else
    if ~isfield(model.ec,'eccodes')
        model.ec.eccodes(1:numel(model.ec.rxns),1) = {''};
    end
    model.ec.eccodes(ecRxns) = eccodes(rxnIdxs(ecRxns));
end
end
