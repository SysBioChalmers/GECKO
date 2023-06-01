classdef GECKOInstaller
% GECKOInstaller
%   Support for installing and uninstalling
%   Run GECKOInstaller.install to install (will set up the path in MATLAB)
%   Run GECKOInstaller.uninstall to clear the path from MATLAB
%   To install, you first need to cd to the GECKO root.

    methods (Static)
        function install
            sourceDir = fileparts(which(mfilename));
            paths = GECKOInstaller.GetFilteredSubPaths(sourceDir, GECKOInstaller.FILE_FILTER);
            addpath(paths);
            savepath;
        end
        function uninstall
            sourceDir = fileparts(which(mfilename));
            paths = GECKOInstaller.GetFilteredSubPaths(sourceDir, GECKOInstaller.FILE_FILTER);
            rmpath(paths);
            savepath;
        end
        function path = getGECKOMainPath()
			path = fileparts(which(mfilename));
			path = strrep(path, '\', '/'); %get rid of backslashes in Windows
			if ~endsWith(path, '/')
				path = strcat(path,'/');
			end
		end

        function newPaths = GetFilteredSubPaths(path_, filter_)
            pathSep = pathsep();
			%Check that there are no separators in the path - that will cause 
            %problems since the separator is used to separate paths in a string
			if contains(path_, pathSep)
				error('The path in which GECKO resides may not contain semicolons for this installation to work!');
			end
            paths = genpath(path_);
            splitPaths = strsplit(paths, pathSep);
            %remove the last, it is empty
            splitPaths = splitPaths(1,1:end-1);
            matches = regexp(splitPaths, filter_, 'match');
            okPaths = cellfun(@isempty, matches);
            pathsLeft = splitPaths(1,okPaths);
            newPaths = strcat(char(join(pathsLeft, pathSep)), pathSep);
        end
    end
    
    properties (Constant)
      FILE_FILTER = '.*\.git|.idea|tutorials.*';
   end
end
