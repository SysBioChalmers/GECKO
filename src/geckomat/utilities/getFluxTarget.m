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
bioRxnIdx    = getIndexes(model, bioRxn,'rxns');
poolIdx      = strcmpi(model.rxns, 'prot_pool_exchange');

% Maximize growth and fix carbon source and suboptimal growth
model = setParam(model, 'obj', bioRxn, 1);
sol   = solveLP(model);
model = setParam(model, 'lb', bioRxn, sol.x(bioRxnIdx) * alpha);

% Minimize prot_pool_exchange and fix
model = setParam(model, 'obj', 'prot_pool_exchange', 1);
sol   = solveLP(model);
model = setParam(model, 'lb', 'prot_pool_exchange', sol.x(poolIdx) * 1.01);

% If minimum fluxes are required get them
minFlux = [];
if nargout == 2
    % Minimize target
    model = setParam(model, 'obj', targetRxn, -1);
    minSol = solveLP(model);
    minFlux = minSol.x;
end

% Maximize target
model = setParam(model, 'obj', targetRxn, 1);
maxSol = solveLP(model);
maxFlux = maxSol.x;
end
