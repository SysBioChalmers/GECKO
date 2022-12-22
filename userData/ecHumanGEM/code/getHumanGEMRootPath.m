function path = getHumanGEMRootPath()
	path = fileparts(which('Human-GEM.mat'));%This file should be on the path
	path = strrep(path, '\', '/'); %get rid of backslashes in Windows
	if ~endsWith(path, '/')
		path = strcat(path,'/');
    end
    %Now remove the model/ at the end
    path = (path(1:strlength(path)-6));
end
