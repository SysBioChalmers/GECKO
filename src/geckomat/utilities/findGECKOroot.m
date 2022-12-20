function [geckoPath, prevDir] = findGECKOroot()
% findGECKOroot
%   Finds the root of the GECKO directory, by searching for the path to
%   GECKO.png. Can also record the current directory, in case a function will
%   use the geckoPath to navigate to a precise folder, and it should return to
%   the previous directory afterwards. 

ST=dbstack('-completenames');
prevDir = pwd();
if length(ST)>1
    geckoPath=ST(2).file; % In case findGECKOroot is run via another function
else
    geckoPath=ST(1).file;
end
rootFound = 0;
while rootFound == 0
    isRoot = exist(fullfile(geckoPath,'GECKO.png'),'file');
    if isRoot == 2
        rootFound = 1;
    else
        ravenPathOld = geckoPath;
        geckoPath = fileparts(geckoPath);
        if strcmp(ravenPathOld,geckoPath)
            error('Cannot find the GECKO root directory. Make sure you have not removed the GECKO.png file from your GECKO installation.')
        end
    end
end