%Abstract Base class for adapters for different species
classdef ModelAdapterManager 
	methods(Static)
        function adapter = getAdapter(adapterPath, addToMatlabPath)
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
                    else
                        printOrange(['WARNING: The adapter will not be on the MATLAB path, since addToMatlabPath is false\n' ...
                                     'and it is not currently on the path. Either set addToMatlabPath to true, fix this\n'...
                                     'manually before calling this function or make sure it is in current directory (not\n'...
                                     'recommended). If the class is not reachable throughout the entire GECKO use there will\n'...
                                     'be errors throughout.\n']);
                    end
                end

            end
            adapter = feval(adapterClassName);
        end
        
        function out = getDefault()
            out = ModelAdapterManager.setGetDefault();
        end
        
        function adapter = setDefault(adapterPath, addToMatlabPath)
            if nargin < 1 || isempty(adapterPath)
                adapter = ModelAdapterManager.setGetDefault(adapterPath);
                return
            end
            if nargin < 2
                addToMatlabPath = true;
            end
            adapter = ModelAdapterManager.setGetDefault(ModelAdapterManager.getAdapter(adapterPath, addToMatlabPath));
        end
        
    end
	methods(Static,Access = private)
        % This is how they recommend defining static variables in Matlab
        function out = setGetDefault(val)
            persistent defaultAdapter; %will be empty initially
            if nargin
                defaultAdapter = val;
            end
            out = defaultAdapter;
        end
    end
end
