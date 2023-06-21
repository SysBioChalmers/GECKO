function ecModel = copyECtoGEM(ecModel, overwrite)
% copyECtoGEM
%   Copies EC numbers from ecModel.ec.eccodes to ecModel.eccodes.
%
% Input:
%   ecModel     an ecModel in GECKO 3 format (with ecModel.ec structure)
%   overwrite   logical that specifies if existing (non-empty) 
%               ecModel.eccodes entries should be overwritten. Empty
%               ecModel.eccodes entries are always overwritten if possible.
%               default false.
%
% Output:
%   ecModel     ecModel with populated model.eccodes
%
% Usage:
%   ecModel = getECfromGEM(ecModel, overwrite)

if nargin < 2 || isempty(overwrite)
    overwrite = false;
end
if ~isfield(ecModel.ec,'eccodes')
    error('ecModel.ec.eccodes does not exist')
end

[a,b] = ismember(ecModel.rxns,ecModel.ec.rxns);
if ~isfield(ecModel,'eccodes')
    ecModel.eccodes = cell(numel(ecModel.rxns),1);
end

if overwrite
    ecModel.eccodes(a) = ecModel.ec.eccodes(b(a));
else % Only replace emptyEcCodes
    emptyEcCodes = cellfun(@isempty, ecModel.eccodes);
    a = a & emptyEcCodes;
    ecModel.eccodes(a) = ecModel.ec.eccodes(b(a));
end
end
