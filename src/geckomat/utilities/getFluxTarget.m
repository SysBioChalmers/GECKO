function [maxFlux, minFlux] = getFluxTarget(model,targetRxn,csRxn,alpha,tolerance,modelAdapter)
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
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   csRxn           rxn ID for the main carbon source uptake reaction.
%   alpha           scalling factor for desired suboptimal growth.
%                   (Optional, default 0.95)
%   tolerance       numerical tolerance for fixing bounds
%                   (Optional, default 1e-4)
%   modelAdapter    a loaded model adapter. (Optional, will otherwise use
%                   the default model adapter)
%
% Output:
%   minFlux         vector of minimum flux rates at minimum target production,
%                   corresponding to model.rxns
%   maxFlux         vector of maximum flux rates at maximum target production,
%                   corresponding to model.rxns
%
% Usage:
%   [maxFlux, minFlux] = simulateGrowth(ecModel,targetRxn,csRxn,alpha,tol)

if nargin < 6 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 5 || isempty(tolerance)
    tolerance = 1e-4;
end

if nargin < 4 || isempty(alpha)
    alpha = 0.95;
end

% Get relevant rxn indexes
targetRxnIdx = getIndexes(model, targetRxn,'rxns');
csRxnIdx = getIndexes(model, csRxn,'rxns');
bioRxnIdx = getIndexes(model, params.bioRxn,'rxns');

% Maximize growth and fix carbon source and suboptimal growth 
model = setParam(model, 'obj', params.bioRxn, 1);
sol = solveLP(model);
model = setParam(model, 'var', csRxn, sol.x(csRxnIdx), tolerance);
model = setParam(model, 'lb', params.bioRxn, sol.x(bioRxnIdx) * (1-tolerance) * alpha);

% If minimum fluxes are required get them
minFlux = [];
if nargout == 2
    % Minimize target
    model = setParam(model, 'obj', targetRxn, -1);
    sol = solveLP(model);

    % Now fix min value for target and minimize proteins usage
    model.lb(targetRxnIdx) = sol.x(targetRxnIdx) * (1-tolerance);
    model = setParam(model, 'obj', 'prot_pool_exchange', 1);
    minSol = solveLP(model,1);
    minFlux = minSol.x;
end

% Maximize target
model = setParam(model, 'obj', targetRxn, 1);
sol = solveLP(model);

% Now fix max value for target and minimize proteins usage
model.lb(targetRxnIdx) = sol.x(targetRxnIdx) * (1-tolerance);
model = setParam(model, 'obj', 'prot_pool_exchange', 1);
maxSol = solveLP(model);
maxFlux = maxSol.x;

end
