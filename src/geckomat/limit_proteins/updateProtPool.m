function ecModel  = updateProtPool(ecModel, Ptot, modelAdapter)
% updateProtPool
%   Updates the protein pool to compensate for measured proteomics data (in
%   model.ec.concs).
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure)
%   Ptot            total protein content in g/gDCW, overwrites the value
%                   from modelAdapter. For instance, condition-specific
%                   fluxData.Ptot from loadFluxData can be used. If nothing
%                   is provided, the modelAdapter value is used.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
% Output:
%   model           an ecModel where model.ec.concs is populated with
%                   protein concentrations
% Usage:
%   ecModel  = updateProtPool(ecModel, Ptot, modelAdapter)

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin < 2 || isempty(Ptot)
    Ptot = params.Ptot;
    disp(['Total protein content used: ' num2str(Ptot) ' [g protein/gDw]'])
end

% Convert Ptot to mg/gDW if provided in g/gDCW (which is default)
if Ptot < 1
    Ptot = Ptot * 1000;
end

Pmeas = sum(ecModel.ec.concs,'omitnan');
if Pmeas == 0
    error('The model has not yet been constrained with proteomics, as ecModel.ec.concs is empty.')
end

Pnew = (Ptot - Pmeas) * params.f;

if Pnew > 0
    PoolRxnIdx = strcmp(ecModel.rxns,'prot_pool_exchange');
    ecModel.lb(PoolRxnIdx) = -Pnew*params.sigma;
    sol = solveLP(ecModel);
    if isempty(sol.x)
        error(['Estimating the remaining protein pool by subtracting the ' ...
               'sum of measured enzyme concentrations (in ecModel.ec.concs) ' ...
               'from the total protein pool (Ptot) does not yield a functional ' ...
               'model.'])
    end
else
    error('The total measured protein mass exceeds the total protein content.')
end
end
