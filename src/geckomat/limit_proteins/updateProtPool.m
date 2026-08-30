function ecModel  = updateProtPool(ecModel, varargin)
% updateProtPool  Update the protein pool to compensate for proteomics.
%
% Only applicable to ecModels where unmeasured enzymes draw from the
% protein pool while measured enzymes (constrained via model.ec.concs) do
% not: shrinks the protein pool exchange upper bound so it represents only
% the remaining, unmeasured protein content (total protein content Ptot
% minus the measured protein mass). For ecModels generated with GECKO
% 3.2.0 or later, where all enzymes draw from the protein pool regardless
% of whether they are measured, this function raises an error instead;
% use setProtPoolSize (see
% https://github.com/SysBioChalmers/GECKO/issues/375).
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% Ptot : double
%     total protein content in g/gDCW, overwrites the value from
%     modelAdapter. For instance, condition-specific fluxData.Ptot from
%     loadFluxData can be used. If nothing is provided, the modelAdapter
%     value is used.
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% ecModel : struct
%     an ecModel where the protein pool exchange reaction's upper bound
%     has been reduced to reflect only the unmeasured protein content.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     ecModel = updateProtPool(ecModel);
%     ecModel = updateProtPool(ecModel, 'Ptot', 0.5);
%
% See also
% --------
% setProtPoolSize, loadFluxData

p = parseGECKOargs(varargin, { ...
    'Ptot',         []; ...
    'modelAdapter', []});
Ptot         = p.Ptot;
modelAdapter = p.modelAdapter;

% Do not run from GECKO version 3.2.0 onwards. This can be recognized by
% prot_usage reactions that are constrained by proteomics and still draw
% from the protein pool
protRxns = find(startsWith(ecModel.rxns,'usage_prot_'));
poolMetIdx = find(strcmp(ecModel.mets,'prot_pool'));
% Selected proteins with proteomics constraints
constProtRxns = ~(ecModel.ub(protRxns)==1000);
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

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(Ptot)
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
    ecModel.ub(PoolRxnIdx) = Pnew*params.sigma;
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
