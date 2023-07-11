function [minFlux, maxFlux] = getFluxTarget(ecModel,targetRxn,csRxn,alpha,tolerance,modelAdapter)
%getFluxTarget
%
%   Function that performs a series of LP optimizations on an ecModel,
%   by first maximizing biomass, then fixing a suboptimal value and 
%   proceeding to protein pool minimization, last a minimization and
%   maximization of a given production target reaction is performed
%   subject to the optimal production level.
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure).
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   csRxn           rxn ID for the main carbon source uptake reaction.
%   alpha           scalling factor for desired suboptimal growth.
%                   (Optional, default 1)
%   tolerance       numerical tolerance for fixing bounds
%                   (Optional, default 1E-6)
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
%   [minFlux, maxFlux] = simulateGrowth(ecModel,targetRxn,csRxn,alpha,tol)

if nargin < 6 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 5 || isempty(tolerance)
    tolerance = 1e-6;
end

if nargin < 4 || isempty(alpha)
    alpha = 1 - tolerance;
end

% Maximize growth
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel);

% Fix carbon source uptake
csRxnIdx = strcmpi(ecModel.rxns,csRxn);
ecModel = setParam(ecModel, 'var', csRxn, sol.x(csRxnIdx), tolerance);

% Fix growth suboptimal
bioRxnIdx = strcmpi(ecModel.rxns, params.bioRxn);
ecModel = setParam(ecModel, 'lb', params.bioRxn, sol.x(bioRxnIdx) * (1-tolerance) * alpha);

% Minimize prot_pool_exchange and fix
poolIdx = strcmpi(ecModel.rxns, 'prot_pool_exchange');
ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
sol = solveLP(ecModel);
ecModel = setParam(ecModel, 'lb', 'prot_pool_exchange', sol.x(poolIdx) * 1.01);

% Minimize target
ecModel = setParam(ecModel, 'obj', targetRxn, -1);
sol = solveLP(ecModel);
minFlux = sol.x;

% Maximize target
ecModel = setParam(ecModel, 'obj', targetRxn, 1);
sol = solveLP(ecModel);
maxFlux = sol.x;

end
