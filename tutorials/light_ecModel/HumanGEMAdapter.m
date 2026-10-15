% HumanGEMAdapter  Model adapter for the light_ecModel tutorial.
%
% The parameters are in model_adapter.toml in this folder, which geckopy
% reads as well. This class only adds getHumanGEMRootPath, so it shows how
% a project keeps custom code next to a TOML adapter: subclass TomlAdapter
% and pass the project folder to its constructor.
%
% See also
% --------
% TomlAdapter, loadAdapterParams
classdef HumanGEMAdapter < TomlAdapter
    methods
        function obj = HumanGEMAdapter()
            obj@TomlAdapter(fullfile(findGECKOroot,'tutorials','light_ecModel'));
        end
    end

    methods(Static)
        function path = getHumanGEMRootPath()
            path = fileparts(which('Human-GEM.mat'));% This file should be on the path
            path = strrep(path, '\', '/'); % get rid of backslashes in Windows
            if ~endsWith(path, '/')
                path = strcat(path,'/');
            end
            % Now remove the model/ at the end
            path = (path(1:strlength(path)-6));
        end
	end
end
