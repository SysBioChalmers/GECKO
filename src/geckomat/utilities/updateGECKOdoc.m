function updateGECKOdoc()
% updateGeckoDoc
%	Updates HTML documentation files for all geckomat GECKO functions. Uses the
%   m2html functions that are provided with RAVEN.
%
% Usage: updateGECKOdoc()

%Make sure that RAVEN-provided m2html is used
ravenDir=findRAVENroot();
path(fullfile(ravenDir,'software','m2html'),path);

%Get the GECKO path
geckoDir=findGECKOroot();

%Remove old doc directory
rmdir(fullfile(geckoDir,'doc'),'s');

%Get a non-redundant list of GECKO subdirectories containing MATLAB
%functions. Absolute paths are not compatible with M2HTML, so convert them
%to the relative paths instead.
geckoDirs=dir(fullfile(geckoDir,'src','geckomat','**/*.m'));
geckoDirs=unique({geckoDirs.folder})';

%Make relative path
relStart = numel(geckoDir)+2;
for i=1:numel(geckoDirs)
    geckoDirs{i,1} = geckoDirs{i,1}(relStart:end);
end

%Save the current working directory and go to RAVEN root directory
originalDir=pwd;
cd(geckoDir);
%Generate HTML documentation files for RAVEN MATLAB functions
m2html('mFiles',geckoDirs,'htmldir','doc');
%Go back to the original working directory
cd(originalDir);
end
