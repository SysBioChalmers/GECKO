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

            GECKOInstaller.checkRAVENversion('2.9.1'); % Minimum RAVEN version

            % Check unique function names
            if ~exist("checkFunctionUniqueness.m")
                error(['Cannot find RAVEN Toolbox in the MATLAB path. Make ' ...
                    'sure you have installed RAVEN in accordance to the ' ...
                    'following instructions, including running ''checkInstallation()'': ' ...
                    'https://github.com/SysBioChalmers/RAVEN/wiki/Installation'])
            else
                status=checkFunctionUniqueness(paths);
                if ~status
                    error(['You might have multiple GECKO installations in your ' ...
                        'MATLAB path. Rerun GECKOInstaller.install after ' ...
                        'resolving the conflicting functions.'])
                end
            end
            addpath(paths);
            savepath;
            GECKOInstaller.checkGECKOversion;
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
            newPaths = char(join(pathsLeft, pathSep));
        end

        function checkRAVENversion(minmVer)
            wrongVersion = false;
            try
                [currVer, installType] = checkInstallation('versionOnly');
                if strcmp(currVer,'develop')
                    printOrange('WARNING: Cannot determine your RAVEN version as it is in a development branch.\n');
                else
                    wrongVersion = GECKOInstaller.isOlderVersion(currVer, minmVer);
                end
            catch
                warning(['Cannot find RAVEN Toolbox in the MATLAB path, or the version ' ...
                    'is too old for this GECKO version (RAVEN ' minmVer ' is required). ' ...
                    'Make sure you have installed RAVEN following the instructions available '...
                    '<a href="https://github.com/SysBioChalmers/RAVEN/wiki/Installation#installation-instructions">here</a>, '...
                    'including running ''checkInstallation()''.'])
            end
            if wrongVersion
                switch installType
                    case 0
                        installType = 'advanced';
                    case 1
                        installType = 'easy';
                    case 2
                        installType = 'medium';
                end
                error(['Installed RAVEN version is %s, while the minimum is %s. '...
                    'Upgrade RAVEN by following the instructions available ' ...
                    '<a href="https://github.com/SysBioChalmers/RAVEN/wiki/Installation#upgrade-raven-after-' ...
                    installType '-installation">here</a>. Do not attempt to run GECKO before upgrading.'], ...
                    currVer,minmVer)
            end
        end

        function checkGECKOversion
            sourceDir = fileparts(which(mfilename));
            hasGit=isfolder(fullfile(sourceDir,'.git'));
            
            if exist(fullfile(sourceDir,'version.txt'), 'file') == 2
                currVer = fgetl(fopen(fullfile(sourceDir,'version.txt')));
                fclose('all');
                fprintf('GECKO version %s installed',currVer)
                try
                    newVer=strtrim(webread('https://raw.githubusercontent.com/SysBioChalmers/GECKO/main/version.txt'));
                    if GECKOInstaller.isOlderVersion(currVer, newVer)
                        fprintf(', newer version %s is available',newVer)
                        if ~hasGit
                            fprintf('\nRun git pull in your favourite git client to update GECKO\n');
                        else
                            fprintf('\nInstructions on how to upgrade <a href="https://github.com/SysBioChalmers/GECKO/wiki/Installation-and-upgrade#installation">here</a>\n');
                        end
                    else
                        fprintf('\n');
                    end
                catch
                    fprintf('\n');
                end
            else
                fprintf('GECKO installed, unknown version (cannot find version.txt).\n')
            end
        end
    end

    methods (Static, Hidden)
        function older = isOlderVersion(ver, refVer)
            % isOlderVersion  True if version string ver precedes refVer.
            % Compares major.minor.patch numerically, from major down to
            % the first number that differs. A pre-release suffix on the
            % patch number (e.g. 4.0.0b1) ranks before the release itself
            % (4.0.0); two pre-releases of one release compare as equal.
            [verNum, verPre] = GECKOInstaller.parseVersion(ver);
            [refNum, refPre] = GECKOInstaller.parseVersion(refVer);
            firstDiff = find(verNum ~= refNum, 1);
            if isempty(firstDiff)
                older = verPre && ~refPre;
            else
                older = verNum(firstDiff) < refNum(firstDiff);
            end
        end

        function [num, isPre] = parseVersion(ver)
            % parseVersion  Leading numbers of major.minor.patch, and
            % whether the patch number carries a pre-release suffix.
            parts = strsplit(strtrim(char(ver)), '.');
            parts(end+1:3) = {'0'};
            num = zeros(1,3);
            for i = 1:3
                lead = regexp(parts{i}, '^\d+', 'match', 'once');
                if ~isempty(lead)
                    num(i) = str2double(lead);
                end
            end
            isPre = ~isempty(regexp(parts{3}, '^\d*[^\d]', 'once'));
        end
    end

    properties (Constant)
        FILE_FILTER = '.*\.git|\.idea|tutorials.*|\.github|_MACOSX|doc';
    end
end
