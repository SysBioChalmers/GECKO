function model = constrainFluxData(model, fluxData, condition, maxMinGrowth, looseStrictFlux, modelAdapter)
% constrainFluxData
%   Constrains fluxes to the data that is provided in the fluxData
%   structure, which itself is read by loadFluxData from data/fluxData.tsv.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   fluxData        structure with flux data
%                   conds       sampling condition
%                   Ptot        total protein (g/gDCW)
%                   grRate      growth rate (1/h)
%                   exchFluxes  exchange fluxes (mmol/gDCWh)
%                   exchMets    exchanged metabolites, matching exchFluxes
%                   exchRxnIDs  exchange reaction IDs, matching exchMets
%   condition       either index number or name of the sample condition in
%                   fluxData.conds (Optional, default = 1)
%   maxMinGrowth    'max' if the provided growth rate should be set as
%                   maximum growth rate (= upper bound), or 'min' if it
%                   should be set as minimum growth rate (= lower bound).
%                   The latter option is suitable if minimization of
%                   prot_pool_exchange is used as objective function. (Opt,
%                   default = 'max')
%   looseStrictFlux how strictly constrained the exchange fluxes should be,
%                   optional, default = 'loose'
%                   'loose' if the exchange fluxes should be constraint
%                           only by the "outer bounds". If exchFluxes(i)
%                           > 0, LB = 0 and UB = exchFluxes(i). If
%                           exchFluxes(i) < 0, LB = exchFluxes(i) and
%                           UB = 0
%                   0-100   LB and UB constraints are set with a specified
%                           percentage of variance around exchFluxes. If 10
%                           is specified, LB = exchFluxes*0.95 and UB =
%                           exchFluxes*1.05. This allows for 10% variance
%                           around the exchFluxes values, but strictly
%                           forces a flux through the exchRxns.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use
%                   the default model adapter)
%
% Output:
%   model           an ecModel where fluxes are constraint
%
% Note: If a provided constraint is either -1000 or 1000, then the function
% will update the reaction lower and upper bound to either allow uptake or
% excretion, irrespective of what option is given as the looseStrictFlux
% parameter.
%
% Usage:
%   model = constrainFluxData(model, fluxData, condition, maxMinGrowth, looseStrictFlux, modelAdapter)

if nargin < 6 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 5 || isempty(looseStrictFlux)
    looseStrictFlux = 'loose';
end

if nargin < 4 || isempty(maxMinGrowth)
    maxMinGrowth = 'max';
end

if nargin < 2 || isempty(fluxData)
    fluxData = loadFluxData(fullfile(params.path,'data','fluxData.tsv'),modelAdapter);
end

if nargin < 3 || isempty(condition)
    condition = 1;
elseif ~isnumeric(condition)
    idx = find(strcmp(fluxData.conds,condition));
    if isempty(condition)
        error(['Condition ' condition ' cannot be found in fluxData'])
    else
        condition = idx;
    end
end

% Select the exchange fluxes of specified condition
fluxData.exchFluxes = fluxData.exchFluxes(condition,:);

%Set original c-source to zero
model = setParam(model,'eq',params.c_source,0);
%Set growth
switch maxMinGrowth
    case 'max'
        model = setParam(model,'lb',params.bioRxn,0);
        model = setParam(model,'ub',params.bioRxn,fluxData.grRate(condition));
    case 'min'
        model = setParam(model,'lb',params.bioRxn,fluxData.grRate(condition));
        model = setParam(model,'ub',params.bioRxn,1000);
end

negFlux = lt(fluxData.exchFluxes,0); % less than 0
ub = fluxData.exchFluxes(~negFlux);
posFlux = fluxData.exchRxnIDs(~negFlux);
lb = fluxData.exchFluxes(negFlux);
negFlux = fluxData.exchRxnIDs(negFlux);

switch looseStrictFlux
    case 'loose'
        model = setParam(model,'lb',negFlux,lb);
        model = setParam(model,'ub',negFlux,0);
        model = setParam(model,'lb',posFlux,0);
        model = setParam(model,'ub',posFlux,ub);
    otherwise
        model = setParam(model,'var',fluxData.exchRxnIDs,fluxData.exchFluxes,looseStrictFlux);
end
extremeFlux = find(abs(fluxData.exchFluxes)==1000);
if any(extremeFlux)
    exchFlux = fluxData.exchFluxes(extremeFlux);
    if any(exchFlux==-1000)
        model = setParam(model,'lb',fluxData.exchRxnIDs(extremeFlux(exchFlux==-1000)),-1000);
        model = setParam(model,'ub',fluxData.exchRxnIDs(extremeFlux(exchFlux==-1000)),0);
    end
    if any(exchFlux==1000)
        model = setParam(model,'lb',fluxData.exchRxnIDs(extremeFlux(exchFlux==1000)),0);
        model = setParam(model,'ub',fluxData.exchRxnIDs(extremeFlux(exchFlux==1000)),1000);
    end
end
end
