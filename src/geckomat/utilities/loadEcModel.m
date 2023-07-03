function model = loadEcModel(filename, modelAdapter)
% loadEcModel
%   Loads the ecModel that matches the modelAdapter. By default, it loads
%   the models/ecModel.yml in the directory specified as param.path in
%   the modelAdapter. Alternative files in the same folder can be loaded by
%   providing the appropriate filename. If loading models from other
%   locations, one can directly use readYAMLmodel.
%
% Input:
%   filename        name of the ecModel file (Optional, default 'ecModel.yml').
%                   Fill should be located in the models/ subfolder of
%                   param.path as specified in the modelAdapter.
%   modelAdapter    a loaded model adapter, from where the model folder is
%                   read (Optional, will otherwise use the default model adapter).
%
% Output:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%
% Usage:
%   model = loadEcModel(filename, modelAdapter)

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
if nargin < 1 || isempty(filename)
    filename = 'ecModel.yml';
else
end

if ~endsWith(filename,{'yml','yaml'})
    error(['ecModels should preferably be distributed in YAML file format, ', ...
           'as otherwise model content will be lost. If the model is in ', ...
           'SBML format, you can use the usual importModel function to load ', ...
           'this as any other genome-scale model.'])
end

filename = fullfile(params.path,'models',filename); 
if endsWith(filename,{'yml','yaml'})
    model = readYAMLmodel(filename);
elseif endsWith(filename,{'xml','sbml'})
    model = importModel(filename);
end
end
