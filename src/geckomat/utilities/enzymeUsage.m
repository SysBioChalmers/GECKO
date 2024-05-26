function usageData = enzymeUsage(ecModel,fluxes,zero)
% enzymeUsage
%   Gives enzyme usages based on a provided flux distribution, as obtained
%   from a full GECKO model. It can give:
%   1)  absolute usage: the specific enzyme usage in mg/gDCW, which can
%       be given for enzymes with- and without concentration information;
%   2)  capacity usage: the ratio of available enzyme that is used, calcuted
%       by (absUsage/UB) (note that capacity usage is 0 if an enzyme
%       concentration was not constrained in the model);
%   3)  UB: the upper bound of each enzyme exchange reaction, which may not
%       be the same as the enzyme concentration, if it has been
%       flexibilized
%   4)  protID: the protein identifiers for each enzyme (if the model has an
%       enzymes field than this order is used, otherwise it is given
%       alphabetically.
%
% Input:
%   ecModel     an ecModel in GECKO 3 format (with ecModel.ec structure)
%   fluxes      vector of fluxes, for instance sol.x
%   zero        logical whether also enzymes with zero absolute usage
%               should be included (Optional, default true)
%
% Output:
%   usageData   structure with enzyme usage data
%               capUsage    vector of enzyme capacity usages
%               absUsage    vector of absolute enzyme usages
%               UB          vector of enzyme exchange reaction upper bounds
%               protID      string array of matching protein IDs
%               fluxes      vector of fluxes, copy of input fluxes
%
% Usage:
%   usageData = enzymeUsage(ecModel,fluxes,zero)

if nargin<3
    zero=true;
end
if ecModel.ec.geckoLight
    error('This function does not work on GECKO light models.')
end
usageData.protID      = ecModel.ec.enzymes;
[~,rxnIdx] = ismember(strcat('usage_prot_',ecModel.ec.enzymes),ecModel.rxns);

usageData.LB          = ecModel.lb(rxnIdx);
usageData.absUsage    = abs(fluxes(rxnIdx));
usageData.capUsage    = abs(usageData.absUsage./usageData.LB);
usageData.fluxes      = fluxes;

if ~zero
    nonzero               = usageData.absUsage<0;
    usageData.absUsage    = usageData.absUsage(nonzero);
    usageData.capUsage    = usageData.capUsage(nonzero);
    usageData.LB          = usageData.LB(nonzero);
    usageData.protID      = usageData.protID(nonzero);
end
end
