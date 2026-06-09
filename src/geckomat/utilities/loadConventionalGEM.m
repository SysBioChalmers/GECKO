function model = loadConventionalGEM(filename, modelAdapter)
% loadConventionalGEM  Load the conventional GEM (non-ecModel).
%
% Loads the conventional GEM (non-ecModel) from the location specified in
% the modelAdapter. By default, it looks in the models/ subdirectory of the
% param.path specified in the modelAdapter. When loading conventional GEMs
% from other locations, one can directly use importModel.
%
% Parameters
% ----------
% filename : char, optional
%     name of the model file, located in the models/ subfolder of
%     param.path as specified in the modelAdapter (default: the value
%     specified as param.convGEM in the modelAdapter).
% modelAdapter : ModelAdapter, optional
%     a loaded model adapter, from where the model folder is read (default:
%     the current default model adapter).
%
% Returns
% -------
% model : struct
%     model in RAVEN format.
%
% Examples
% --------
%     model = loadConventionalGEM(filename, modelAdapter);
%
% See also
% --------
% loadEcModel


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
