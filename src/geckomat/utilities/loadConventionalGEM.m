function model = loadConventionalGEM(varargin)
% loadConventionalGEM  Load the conventional GEM (non-ecModel).
%
% Loads the conventional GEM (non-ecModel) from the location specified in
% the modelAdapter. By default, it looks in the models/ subdirectory of the
% param.path specified in the modelAdapter. When loading conventional GEMs
% from other locations, one can directly use importModel.
%
% Name-Value Arguments
% --------------------
% filename : char
%     name of the model file, without extension, located in the models/
%     subfolder of param.path as specified in the modelAdapter; it is
%     loaded as SBML (.xml) (default: the value specified as param.convGEM
%     in the modelAdapter, which may instead point to a YAML file).
% modelAdapter : ModelAdapter
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
%     % optional arguments may be given positionally or as name-value pairs:
%     model = loadConventionalGEM(filename, modelAdapter);
%     model = loadConventionalGEM('filename', filename);
%
% See also
% --------
% loadEcModel

p = parseGECKOargs(varargin, { ...
    'filename',     []; ...
    'modelAdapter', []});
filename     = p.filename;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
if isempty(filename)
    filename = params.convGEM;
else
    filename = fullfile(params.path,'models',[filename, '.xml'])
end

if endsWith(filename,'.xml')
    model = importModel(filename);
elseif endsWith(filename,'.yml')
    model = readYAMLmodel(filename);
end
