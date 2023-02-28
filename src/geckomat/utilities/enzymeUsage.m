function [absUsage,capUsage,UB,protID] = enzymeUsage(ecModel,fluxes,zero)
% enzymeUsage
%   Gives enzyme usages based on a provided flux distribution. It can give:
%   1)  absolute usage: the specific enzyme usage in ug/gDCW/h, which can
%       be given for enzymes with- and without concentration information;
%   2)  capacity usage: the ratio of available enzyme that is used, calcuted
%       by (absUsage/UB) (note that capacity usage is 0 if an enzyme
%       concentration was not constrained in the model);
%   3)  UB: the upper bound of each enzyme exchange reaction, which may not
%       be the same as the enzyme concentration, if it has been
%       flexibilized
%   4)  protId: the protein identifiers for each enzyme (if the model has an
%       enzymes field than this order is used, otherwise it is given
%       alphabetically.
%
%  Input:
%   ecModel         a GECKO3 ecModel
%   fluxes          vector of fluxes, for instance sol.x
%   zero            logical whether also enzymes with zero absolute usage
%                   should be included (Optional, default true)
%
%  Output:
%   capUsage        vector of enzyme capacity usages
%   absUsage        vector of absolute enzyme usages
%   UB              vector of enzyme exchange reaction upper bounds
%   protID          string array of protein IDs matching the other output
%
% Usage: [absUsage,capUsage,UB,protID] = enzymeUsage(ecModel,fluxes,zero)

if nargin<3
    zero=true;
end

protID      = ecModel.ec.enzymes;
[~,rxnIdx] = ismember(strcat('usage_prot_',ecModel.ec.enzymes),ecModel.rxns);
%[~,rxnIdx]  = ismember(ecModel.rxns,strcat('usage_prot_',ecModel.ec.enzymes));

UB          = ecModel.ub(rxnIdx);
absUsage    = fluxes(rxnIdx);
capUsage    = absUsage./UB;

if ~zero
    nonzero     = absUsage>0;
    absUsage    = absUsage(nonzero);
    capUsage    = capUsage(nonzero);
    UB          = UB(nonzero);
    protId      = protId(nonzero);
end
end
