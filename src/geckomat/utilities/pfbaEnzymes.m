function solution = pfbaEnzymes(model, varargin)
% pfbaEnzymes  Parsimonious FBA minimising total enzyme usage.
%
% Among all flux distributions at (or within a fraction of) the optimal
% value of the current objective, finds the one that minimises total
% enzyme usage: the smallest sum of usage_prot_* fluxes. This is the
% enzyme-aware analogue of classical pFBA, which minimises total flux
% instead --- for ecModels this usually picks the most parsimonious
% *proteome*, not just the most parsimonious fluxes.
%
% Parameters
% ----------
% model : struct
%     a full ecModel in GECKO 3 format (with ecModel.ec structure and
%     usage_prot_<id> reactions). Not supported for gecko-light models.
%
% Name-Value Arguments
% --------------------
% rxnId : char
%     reaction id to use as the objective instead of model.c as given
%     (default: use model.c as given).
% fractionOfOptimum : double
%     fraction of the optimal objective value to fix before minimising
%     enzyme usage (default 1.0).
%
% Returns
% -------
% solution : struct
%     - x : flux distribution (indexed by model.rxns) at the
%       enzyme-usage-minimising solution.
%     - enzymeUsage : the minimised sum of usage_prot_* fluxes.
%     - objectiveValue : the fixed objective's own value at this solution.
%     - stat : solveLP exit flag.
%
% Examples
% --------
%     solution = pfbaEnzymes(ecModel);
%     solution = pfbaEnzymes(ecModel, 'fractionOfOptimum', 0.9);
%
% Notes
% -----
% Ported from geckopy's pfba_enzymes (utilities/pfba_enzymes.py), itself
% ported from the legacy geckopy package (Carrasco et al., 2023,
% https://doi.org/10.1128/spectrum.01705-23), geckopy/flux_analysis.py:342-386
% (pfba_protein).
%
% Fixes the objective as a constraint the same way solveLP's own minFlux=1
% option does internally: a "fake metabolite" produced by each reaction
% with the stoichiometry it has in the objective, constrained to at least
% the target value (with the same 1e-6 relative safety margin solveLP
% itself applies, to avoid spurious infeasibility from solver tolerance).
%
% usage_prot_<id> reactions are forward-only (lb=0) in the standard GECKO 3
% layout. If any has lb<0 (reverse flux enabled), minimising its raw flux
% would not minimise |flux| for that reaction, unlike geckopy's cobra-based
% forward/reverse variable split --- this function errors in that case
% rather than silently building an incorrect objective.
%
% See also
% --------
% getEnzymeBottlenecks

p = parseGECKOargs(varargin, {'rxnId', []; 'fractionOfOptimum', []});
rxnId             = p.rxnId;
fractionOfOptimum = p.fractionOfOptimum;
if isempty(fractionOfOptimum)
    fractionOfOptimum = 1.0;
end

if isfield(model.ec, 'geckoLight') && model.ec.geckoLight
    error(['pfbaEnzymes requires a full ecModel, not gecko-light: light models have no ' ...
        'usage_prot_<id> reactions.'])
end

if ~isempty(rxnId)
    rxnIdx = find(strcmp(model.rxns, rxnId));
    if isempty(rxnIdx)
        error('pfbaEnzymes: reaction %s not found in model.rxns.', rxnId)
    end
    model.c = zeros(numel(model.rxns), 1);
    model.c(rxnIdx) = 1;
end
if ~any(model.c)
    error('pfbaEnzymes: model.c is empty; nothing to fix before minimising enzyme usage.')
end

usageIdx = find(startsWith(model.rxns, 'usage_prot_'));
if isempty(usageIdx)
    error('pfbaEnzymes: no usage_prot_<id> reactions found. This is done by makeEcModel.')
end
if any(model.lb(usageIdx) < 0)
    error(['pfbaEnzymes does not support usage_prot_<id> reactions with a negative lower ' ...
        'bound (reverse flux enabled): minimising raw flux would not minimise |flux| for ' ...
        'those reactions.'])
end

solOrig = solveLP(model);
if solOrig.stat ~= 1
    error('pfbaEnzymes: solver status %d on the original objective; cannot fix it as a constraint.', solOrig.stat)
end

% Fix the objective as a constraint, same "fake metabolite" technique
% solveLP's own minFlux=1 option uses (see solveLP.m).
target = fractionOfOptimum * solOrig.f;
if target > 0
    target = target * 0.999999;
else
    target = target * 1.000001;
end
fixModel = model;
fixModel.S(end+1, :) = model.c';
fixModel.mets{end+1, 1} = 'pfbaEnzymes_objective';
if size(fixModel.b, 2) == 1
    fixModel.b = [fixModel.b fixModel.b];
end
fixModel.b(end+1, :) = [target inf];

fixModel.c = zeros(numel(fixModel.rxns), 1);
fixModel.c(usageIdx) = -1;

solution = solveLP(fixModel);
if solution.stat ~= 1
    error('pfbaEnzymes: minimising enzyme usage was infeasible at fractionOfOptimum=%g.', fractionOfOptimum)
end

solution.enzymeUsage    = -solution.f;
solution.objectiveValue = model.c' * solution.x;
end
