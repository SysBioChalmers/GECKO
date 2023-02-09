function model = constrainPool(model, Ptot, f, sigma, modelAdapter)
% constrainPool
%   Draw all enzyme usages from the total protein pool, discarding any
%   constraints based on specific protein concentrations.
%
% Input:
%   model           an ec-model
%   Ptot            Total cellular protein content in g/gDCW. If not
%                   specified, the value will be read from the model
%                   adapter. If not specified in model adapter, 0.5 g/gDCW
%                   is assumed.
%   f               Estimated fraction of enzymes in the model. If not
%                   specified, the value will be read from the model
%                   adapter. If not specified in model adapter, 0.5 is
%                   assumed.
%   sigma           Estimated saturation factor. If not specified, the
%                   value will be read from the model adapter. If not
%                   specified in model adapter, 0.5 is assumed.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
% Output:
%   model           an ec-model constraint by total protein

%Get indices of usage reactions 
usageRxns = strcat('usage_prot_',model.ec.enzymes);
[~, usageRxnsIdx] = ismember(usageRxns, model.rxns);

if any(usageRxnsIdx == 0)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end
%Get index of protein pool exchange rxn
protPoolIdx = find(ismember(model.mets,'prot_pool'));
if ~any(protPoolIdx)
    error('Cannot find protein pool pseudometabolite.')
end

%Set all reactions to draw from prot_pool
model.S(protPoolIdx, usageRxnsIdx) = -1;
model.ub(usageRxnsIdx) = Inf;

if nargin<5
    modelAdapter = [];
end
if nargin<4
    sigma = [];
end
if nargin<3
    f = [];
end
if nargin<2
    Ptot = [];
end

model = setProtPoolSize(model, Ptot, f, sigma, modelAdapter);
end
