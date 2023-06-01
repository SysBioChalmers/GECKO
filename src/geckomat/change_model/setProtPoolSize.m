function model = setProtPoolSize(model, Ptot, f, sigma, modelAdapter)
% setProtPoolSize
%   Sets the limit of the total protein usage in the model.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
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
%   model           ecModel with protein pool constraint set
%
% Usage: model = setProtPoolSize(model, Ptot, f, sigma, modelAdapter)

if nargin < 5 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter) && (nargin<2 || isempty(Ptot) || isempty(f) || isempty(sigma))
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if nargin<4 || isempty(sigma)
    sigma = modelAdapter.getParameters().sigma;
end
if nargin<3 || isempty(f)
    f = modelAdapter.getParameters().f;
end
if nargin<2 || isempty(Ptot)
    Ptot = modelAdapter.getParameters().Ptot;
end

model.lb(strcmp(model.rxns, 'prot_pool_exchange')) = -(Ptot*f*sigma*1000);
end

