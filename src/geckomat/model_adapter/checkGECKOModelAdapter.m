function params = checkGECKOModelAdapter(GECKOModelAdapter)
% checkGECKOModelAdapter
%   Checks whether GECKOModelAdapter is set as global variable (if not, the
%   user should run setGECKOModelAdapter). If an output is specified, this
%   function gives the model-specific parameters as specified in the
%   ModelAdapter file under getParameters. This function should be run at
%   the start of each function that requires information from ModelAdapter.
%
% Input:
%   GECKOModelAdapter   the GECKOModelAdapter parameter that should have
%                       been set by setGECKOModelAdapter.
% 
% Usage: params = checkGECKOModelAdapter(GECKOModelAdapter)
 
if isempty(GECKOModelAdapter)
    error(['GECKO ModelAdapter is not set. Prepare a model-specific '...
           'ModelAdapter (see userData/ecYeastGEM/ModelAdapter.m for an '...
           'example) and set this as the ModelAdapter with '...
           'setGECKOmodelAdapter'])
end
if nargout>0
    params=GECKOModelAdapter.getParameters();
end
