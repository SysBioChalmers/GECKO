function model = constrainEnzConcs(model)
% constrainEnzConcs
%   Constrain enzyme usages by their concentration as provided in
%   model.ec.concs. For enzymes with non-NaN entries in model.ec.concs,
%   their enzyme usage reaction will no longer draw from the protein pool,
%   but is rather constraint by the measured protein abundance.
%
% Input:
%   model   an ecModel in GECKO 3 format (with ecModel.ec structure) with enzyme
%           concentrations in model.ec.concs
%
% Output:
%   model   an ecModel constraint with available enzyme concentrations
%
% Note: to populate model.ec.concs you should run getProteomics.
%
% Usage:
%   model = constrainEnzConcs(model)

%Enzyme with NaN entry in model.ec.concs => draw from prot_pool
%Enzyme with numeric entry in model.ec.concs => exchange reaction with
%enzyme level as UB

%Get indices of usage reactions 
usageRxns = strcat('usage_prot_',model.ec.enzymes);
[~, usageRxnsIdx] = ismember(usageRxns, model.rxns);

if any(usageRxnsIdx == 0)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end
%Get index of protein pool metabolite
protPoolIdx = find(ismember(model.mets,'prot_pool'));
if ~any(protPoolIdx)
    error('Cannot find protein pool pseudometabolite.')
end

%Protein that should be constraint by UB
protCons = ~isnan(model.ec.concs);

%Set all reactions to draw from prot_pool
model.S(protPoolIdx, usageRxnsIdx(~protCons)) = 1;
model.lb(usageRxnsIdx(protCons)) = -1000;

%If non-NaN in model.ec.concs, then constrain by UB
if any(protCons)
    model.S(protPoolIdx, usageRxnsIdx(protCons)) = 0;
    model.lb(usageRxnsIdx(protCons)) = -model.ec.concs(protCons);
end
end
