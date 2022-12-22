function setGECKOModelAdapter(path)
% setGECKOModelAdapter
%   Set the global GECKOModelAdapter to the provided model-specific
%   ModelAdapter. This is required to be set before working on any model
%   with GECKO, and sets model-specific parameters and functions.
%
% Input:
%   path    full file path where the ModelAdapter can be found. Example:
%           'userData/ecYeastGEM/YeastGEMAdapter.m'. If the path starts
%           with 'userData/', it is assumed that this folder is localized
%           in the GECKO directory. Note that the ModelAdapter file cannot
%           be named ´ModelAdapter.m`.
% 
% Usage: setGECKOModelAdapter(path)

if ~exist(path,'file')
    error([path ' file not found.'])
elseif startsWith(path,'userData/') % Likely in GECKO folder
    geckoPath = findGECKOroot;
    path = fullfile(geckoPath,path);
end

global GECKOModelAdapter;
evalc('run(path);');
GECKOModelAdapter = ans;
disp([path ' is set as GECKOModelAdapter'])
end
