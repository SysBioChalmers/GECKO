% ModelAdapterManager  Resolve and cache the active ModelAdapter.
%
% Static helper class that loads a model adapter from a file path and
% stores a process-wide default adapter, so that GECKO functions can
% retrieve the active adapter without it being passed explicitly.
%
% Notes
% -----
% All methods are static. The default adapter is held in a persistent
% variable and shared across the whole GECKO session.
%
% See also
% --------
% ModelAdapter
classdef ModelAdapterManager 
	methods(Static)
        function adapter = getAdapter(adapterPath, addToMatlabPath)
            if nargin < 2
                addToMatlabPath = true;
            end
            ModelAdapterManager.showGECKO4Notice();

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

        function showGECKO4Notice()
            % Remove this method (and its call in getAdapter) in a future
            % GECKO 4.x release, once users have had time to notice the
            % backward-incompatible changes from GECKO 3.
            persistent shown;
            if isempty(shown)
                printOrange(['NOTE: This is GECKO 4, which is not fully backward-compatible\n' ...
                    'with GECKO 3. See https://gecko-docs.readthedocs.io/en/latest/gecko3-to-gecko4.html\n' ...
                    'for what changed and how to stay on GECKO 3 (v3.2.5) if needed.\n']);
                shown = true;
            end
        end
    end
end
