function usageData = getEnzymeUsage(ecModel,fluxes,varargin)
% getEnzymeUsage  Compute enzyme usages from a provided flux distribution.
%
% Gives enzyme usages based on a provided flux distribution, as obtained
% from a full GECKO model. It can give:
%
% 1. absolute usage: the specific enzyme usage in mg/gDCW, which can be
%    given for enzymes with- and without concentration information;
% 2. capacity usage: the ratio of available enzyme that is used, calculated
%    by (absUsage/UB) (note that capacity usage is 0 if an enzyme
%    concentration was not constrained in the model);
% 3. UB: the upper bound of each enzyme exchange reaction, which may not be
%    the same as the enzyme concentration, if it has been flexibilized;
% 4. protID: the protein identifiers for each enzyme, in the same order as
%    ecModel.ec.enzymes.
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% fluxes : double
%     vector of fluxes, for instance sol.x.
%
% Name-Value Arguments
% --------------------
% zero : logical
%     whether also enzymes with zero absolute usage should be included
%     (default true).
%
% Returns
% -------
% usageData : struct
%     structure with enzyme usage data.
%
% Notes
% -----
% The usageData structure has the following fields:
%
% - capUsage : vector of enzyme capacity usages.
% - absUsage : vector of absolute enzyme usages.
% - UB : vector of enzyme exchange reaction upper bounds.
% - protID : string array of matching protein IDs.
% - fluxes : vector of fluxes, copy of input fluxes.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     usageData = getEnzymeUsage(ecModel,fluxes,zero);
%     usageData = getEnzymeUsage(ecModel,fluxes,'zero',false);
%
% See also
% --------
% reportEnzymeUsage

p = parseGECKOargs(varargin, { ...
    'zero', true});
zero = p.zero;
if ecModel.ec.geckoLight
    error('This function does not work on GECKO light models.')
end
usageData.protID      = ecModel.ec.enzymes;
[~,rxnIdx] = ismember(strcat('usage_prot_',ecModel.ec.enzymes),ecModel.rxns);

usageData.UB          = ecModel.ub(rxnIdx);
usageData.absUsage    = abs(fluxes(rxnIdx));
usageData.capUsage    = abs(usageData.absUsage./usageData.UB);
usageData.fluxes      = fluxes;

if ~zero
    nonzero               = usageData.absUsage>0;
    usageData.absUsage    = usageData.absUsage(nonzero);
    usageData.capUsage    = usageData.capUsage(nonzero);
    usageData.UB          = usageData.UB(nonzero);
    usageData.protID      = usageData.protID(nonzero);
end
end
