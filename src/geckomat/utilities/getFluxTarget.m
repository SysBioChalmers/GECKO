function [maxFlux, minFlux] = getFluxTarget(model,bioRxn,targetRxn,alpha)
% getFluxTarget
%
%   Function that performs a series of LP optimizations on an ecModel,
%   by first maximizing biomass, then fixing a suboptimal value and 
%   proceeding to protein pool minimization, last a minimization and
%   maximization of a given production target reaction is performed
%   subject to the optimal production level.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure).
%   bioRxn          rxn ID for the biomass reaction.
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   alpha           scalling factor for desired suboptimal growth.
%                   (Optional, default 0.95)
%
% Output:
%   minFlux         vector of minimum flux rates at minimum target production,
%                   corresponding to model.rxns
%   maxFlux         vector of maximum flux rates at maximum target production,
%                   corresponding to model.rxns
%
% Usage:
%   [maxFlux, minFlux] = getFluxTarget(ecModel,bioRxn,targetRxn,alpha)

if nargin < 4 || isempty(alpha)
    alpha = 0.95;
end

% Fix a suboptimal alpha if equal to 1
if alpha == 1
	alpha = 0.99;
end

% Get relevant rxn indexes
bioRxnIdx    = getIndexes(model, bioRxn, 'rxns');
targetRxnIdx = getIndexes(model, targetRxn, 'rxns');

% Maximize growth and fix carbon source and suboptimal growth
model = setParam(model, 'obj', bioRxn, 1);
sol   = solveLP(model);
model = setParam(model, 'lb', bioRxn, sol.x(bioRxnIdx) * alpha);

% If minimum fluxes are required get them
minFlux = [];
if nargout == 2
    % Minimize target
    model   = setParam(model, 'obj', targetRxn, -1);
    minSol  = solveLP(model);
    % Fix ub to minimum attainable product
    model.ub(targetRxnIdx) = minSol.x(targetRxnIdx) * 1.1;
    % Minimize prot_pool_exchange
    model = setParam(model, 'obj', 'prot_pool_exchange', 1);
    minUsage = solveLP(model);
    minFlux  = minUsage.x;
end

% Maximize target
model = setParam(model, 'obj', targetRxn, 1);
maxSol = solveLP(model);
% Fix lb to 90% of the maximum attainable product. Based on
% https://doi.org/10.1128/AEM.00115-10
model.lb(targetRxnIdx) = maxSol.x(targetRxnIdx) * 0.90;
% Minimize prot_pool_exchange
model = setParam(model, 'obj', 'prot_pool_exchange', 1);
maxUsage = solveLP(model);
maxFlux = maxUsage.x;
end
