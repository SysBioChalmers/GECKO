function bottlenecks = getEnzymeBottlenecks(model, varargin)
% getEnzymeBottlenecks  Top-N enzymes limiting the current objective.
%
% Solves the model and ranks every enzyme by the absolute shadow price of
% its prot_<id> mass-balance constraint --- the enzymes whose extra
% availability would help the objective the most.
%
% Parameters
% ----------
% model : struct
%     a full ecModel in GECKO 3 format (with ecModel.ec structure). Not
%     supported for gecko-light models, which have no per-enzyme
%     prot_<id> / usage_prot_<id> machinery.
%
% Name-Value Arguments
% --------------------
% top : double
%     number of enzymes to return (default 10).
%
% Returns
% -------
% bottlenecks : table
%     one row per returned enzyme, sorted by descending absolute shadow
%     price, with columns uniprot, gene, shadowPrice, flux, capUsage
%     (abs(flux)/upperBound, NaN if upperBound is 0) and upperBound.
%
% Examples
% --------
%     bottlenecks = getEnzymeBottlenecks(ecModel);
%     bottlenecks = getEnzymeBottlenecks(ecModel, 'top', 5);
%
% Notes
% -----
% Ported from geckopy's get_enzyme_bottlenecks (utilities/bottlenecks.py),
% itself ported from the legacy geckopy package (Carrasco et al., 2023,
% https://doi.org/10.1128/spectrum.01705-23), geckopy/flux_analysis.py:275-281
% (get_protein_bottlenecks).

p = parseGECKOargs(varargin, {'top', []});
top = p.top;
if isempty(top)
    top = 10;
end

if isfield(model.ec, 'geckoLight') && model.ec.geckoLight
    error(['getEnzymeBottlenecks requires a full ecModel, not gecko-light: ' ...
        'light models have no per-enzyme prot_<id> / usage_prot_<id> machinery.'])
end

sol = solveLP(model);
if sol.stat ~= 1
    error('getEnzymeBottlenecks: solver status %d; cannot rank bottlenecks.', sol.stat)
end

uniprot = model.ec.enzymes(:);
gene    = model.ec.genes(:);

protMets  = strcat('prot_', uniprot);
[found, metIdx] = ismember(protMets, model.mets);
if ~all(found)
    error('Cannot find prot_<id> metabolites for all enzymes. This is done by makeEcModel.')
end
usageRxns = strcat('usage_prot_', uniprot);
[found, rxnIdx] = ismember(usageRxns, model.rxns);
if ~all(found)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end

shadowPrice = sol.sPrice(metIdx);
flux        = sol.x(rxnIdx);
upperBound  = model.ub(rxnIdx);
capUsage    = abs(flux) ./ upperBound;
capUsage(upperBound == 0) = NaN;

bottlenecks = table(uniprot, gene, shadowPrice, flux, capUsage, upperBound);
[~, order] = sort(abs(bottlenecks.shadowPrice), 'descend');
bottlenecks = bottlenecks(order, :);
bottlenecks = bottlenecks(1:min(top, height(bottlenecks)), :);
end
