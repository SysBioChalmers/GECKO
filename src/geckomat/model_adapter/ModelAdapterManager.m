%Abstract Base class for adapters for different species
classdef ModelAdapterManager 
	methods(Static)
        %not really needed, we could access it directly. Not very nice, but practical
        function adapter = getAdapterFromPath(adapterPath, addToMatlabPath)
            if nargin < 2
                addToMatlabPath = true;
            end
            
            potentialAdapterFiles = dir(fullfile(adapterPath,'*.m')); %gets all wav files in struct
            %Check that there is only one .m file
            if length(potentialAdapterFiles) ~= 1
                error('getAdapterFromPath: One, and only one .m file is expected in the folder. This file is expected to contain a ModelAdapter.');
            end
            %Check if the folder is on the path
            s = pathsep;
            pathStr = [s, path, s];
            onPath = contains(pathStr, [s, adapterPath, s], 'IgnoreCase', ispc);
            
            if ~onPath %otherwise do nothing
               if addToMatlabPath
                   addpath(adapterPath);
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
            
            adapterClassName = extractBefore(potentialAdapterFiles.name,length(potentialAdapterFiles.name)-1);
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
                addToMatlabPath = [];
            end
            ModelAdapterManager.setDefaultAdapter(ModelAdapterManager.getAdapterFromPath(adapterPath, addToMatlabPath));
        end
        
    end
	methods(Static,Access = private)
        %This is how they recommend defining static variables in matlab. Strange, but works
        function out = setGetDefaultAdapter(val)
            persistent defaultAdapter; %will be empty initially
            if nargin
                defaultAdapter = val;
            end
            out = defaultAdapter;
        end
    end

    %To have the params public is a bit "ugly", but very practical 
    %if we want to change a parameter
    properties (Access = private) 
        params;
    end
end
