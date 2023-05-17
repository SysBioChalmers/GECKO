%Abstract Base class for adapters for different species
classdef ModelAdapterManager 
	methods(Static)
        function adapter = getAdapterFromPath(adapterPath, addToMatlabPath)
            if nargin < 2
                addToMatlabPath = true;
            end
            
            [adapterFolder, adapterClassName, extension] = fileparts(adapterPath);
            if ~strcmp(extension, '.m')
                error('Please provide the full path to the adapter file, including the file extension.');
            else
                s = pathsep;
                pathStr = [s, path, s];
                onPath = contains(pathStr, [s, adapterFolder, s], 'IgnoreCase', ispc);
                % Check if the folder is on the path
                if ~onPath
                    if addToMatlabPath
                        addpath(adapterFolder);
                        savepath();
                    else
                        warning(['The adapter will not be on the MATLAB path, since ' ...
                                'addToMatlabPath is false and it is not currently on ' ...
                                'the path. Either set addToMatlabPath to true, fix this' ...
                                'manually before calling this function or make sure it ' ...
                                'is in current directory (not recommended). If the ' ...
                                'class is not reachable throughout the entire GECKO use ' ...
                                'there will be an error.']);
                    end
                end

            end
            adapter = feval(adapterClassName);
        end
        
        function setDefaultAdapter(val)
            ModelAdapterManager.setGetDefaultAdapter(val);
        end
        
        function out = getDefaultAdapter()
            out = ModelAdapterManager.setGetDefaultAdapter();
        end
        
        function setDefaultAdapterFromPath(adapterPath, addToMatlabPath)
            if nargin < 2
                addToMatlabPath = true;
            end
            ModelAdapterManager.setDefaultAdapter(ModelAdapterManager.getAdapterFromPath(adapterPath, addToMatlabPath));
        end
        
    end
	methods(Static,Access = private)
        % This is how they recommend defining static variables in Matlab
        function out = setGetDefaultAdapter(val)
            persistent defaultAdapter; %will be empty initially
            if nargin
                defaultAdapter = val;
            end
            out = defaultAdapter;
        end
    end
end
