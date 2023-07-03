function model = loadConventionalGEM(filename, modelAdapter)
% loadConventionalGEM
%   Loads the conventional GEM (non-ecModel) from the location specified in
%   the modelAdapter. By default, it looks in the models/ subdirectory of 
%   the param.path specified in the modelAdapter. When loading
%   conventional GEMs from other locations, one can directly use importModel.
%
% Input:
%   filename        name of the model file, located in the the models/
%                   subfolder of param.path as specified in the
%                   modelAdapter (Optional, will otherwise use the
%                   value specified as param.convGEM in the
%                   modelAdapter)
%   modelAdapter    a loaded model adapter, from where the model folder is
%                   read (Optional, will otherwise use the default model adapter).
%
% Output:
%   model           model in RAVEN format
%
% Usage:
%   model = loadConventionalGEM(filename, modelAdapter)


if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
if nargin < 1 || isempty(filename)
    filename = params.convGEM;
else
    filename = fullfile(params.path,'models',[filename, '.xml'])
end

if endsWith(filename,'.xml')
    model = importModel(filename);
elseif endsWith(filename,'.yml')
    model = readYAMLmodel(filename);
end
