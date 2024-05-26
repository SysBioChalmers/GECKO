function ecModel  = updateProtPool(ecModel, Ptot, modelAdapter)
% updateProtPool
%   Obsolete since GECKO 3.2.0, as all (measured and unmeasured) enzymes
%   are drawn from the protein pool. Instead, use setProtPoolSize. See
%   https://github.com/SysBioChalmers/GECKO/issues/375 for explanation.
%
%   Before GECKO 3.2.0: updates the protein pool to compensate for measured
%   proteomics data (in model.ec.concs), as only the unmeasured enzymes
%   draw from the protein pool.
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure)
%   Ptot            total protein content in g/gDCW, overwrites the value
%                   from modelAdapter. For instance, condition-specific
%                   fluxData.Ptot from loadFluxData can be used. If nothing
%                   is provided, the modelAdapter value is used.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
% Output:
%   model           an ecModel where model.ec.concs is populated with
%                   protein concentrations
% Usage:
%   ecModel  = updateProtPool(ecModel, Ptot, modelAdapter)

% Do not run from GECKO version 3.2.0 onwards. This can be recognized by
% prot_usage reactions that are constrained by proteomics and still draw
% from the protein pool
protRxns = find(startsWith(ecModel.rxns,'usage_prot_'));
poolMetIdx = find(strcmp(ecModel.mets,'prot_pool'));
% Selected proteins with proteomics constraints
constProtRxns = ~(ecModel.lb(protRxns)==-1000);
% Are any still drawing from prot_pool? This is introduced in GECKO 3.2.0.
if any(full(ecModel.S(poolMetIdx,protRxns(constProtRxns))))
    error(['In the provided ecModel, all protein usage reactions, both with ' ...
        'and without concentration constraints, draw from the protein pool. ' ...
        'This was introduced with GECKO 3.2.0. Since then, updateProtPool ' ...
        'has become obsolete, use setProtPoolSize instead to constrain the ' ...
        'total protein pool with  condition-specific total protein content. See '...
        '<a href="https://github.com/SysBioChalmers/GECKO/issues/375">here</a> ' ...
        'for more explanation.'])
end

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
    error('The model has not yet been populated with proteomics, as ecModel.ec.concs is empty.')
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
