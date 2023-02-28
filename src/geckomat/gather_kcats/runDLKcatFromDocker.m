function runDLKcatFromDocker(deleteImage, DLKcatFile, modelAdapter, DLKcatPath)
% runDLKcatFromDocker
%   Runs DLKcat to predict kcat values from a Docker image
%
% Input
%   deleteImage      true or false if delete the docker image tar.gz file.
%                   (Optional, default = false)
%   DLKcatFile      path to the DLKcat.tsv file (including file name), as
%                   written by writeDLKcatFile. Once DLKcat is succesfully
%                   run, the DLKcatFile will be overwritten with the DLKcat
%                   output. Optional, otherwise the file location will be
%                   assumed to be in the model-specific 'data' sub-folder
%                   taken from modelAdapter (e.g.
%                   GECKO/userData/ecYeastGEM/data/DLKcat.tsv)
%   modelAdapter    a loaded model adapter. (Optional, will otherwise use
%                   the default model adapter)
%   DLKcatPath      path where DLKcat is/will be installed. (Optional,
%                   defaults to GECKO/dlkcat)
%
%   NOTE: 1. Requires Docker to be installed, and docker desktop running. Visit "https://www.docker.com"
%         2. Runtime will depend on whether the image is to be downloaded or not.

%Get the GECKO path
geckoPath = findGECKOroot();

if nargin < 4 || isempty(DLKcatPath)
    DLKcatPath = fullfile(geckoPath,'dlkcat');
end

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin < 2 || isempty(DLKcatFile)
    DLKcatFile = fullfile(params.path,'data','DLKcat.tsv');
end

if nargin < 1 || isempty(DLKcatPath)
    deleteImage = false;
end

%% Check and install requirements
% On Mac, docker might not be properly loaded if MATLAB is started via
% launcher and not terminal.
if ismac
    setenv('PATH', strcat('/usr/local/bin', ':', getenv("PATH")));
end

% Check if docker is installed
[checks.docker.status, checks.docker.out] = system('docker --version');
if checks.docker.status ~= 0
    error('Cannot find Docker. Make sure it is installed')
end
webOptions = weboptions('Timeout',30);
% Dowload DLKcat package
if ~exist(fullfile(DLKcatPath,'DLKcat.py'),'file')
    if ~exist(fullfile(DLKcatPath),'dir')
        mkdir(fullfile(DLKcatPath));
    end
    disp('=== Downloading DLKcat ...')
    packageURL = 'https://github.com/SysBioChalmers/GECKO/raw/dlkcatPackage/dlkcat.zip';
    %packageURL = 'https://github.com/SysBioChalmers/GECKO/releases/download/v3.0.0/dlkcat_package.zip';
    websave(fullfile(DLKcatPath,'dlkcat_package.zip'),packageURL,webOptions);
    unzip(fullfile(DLKcatPath,'dlkcat_package.zip'),geckoPath);
    delete(fullfile(DLKcatPath,'dlkcat_package.zip'));
end

currPath = pwd();
cd(DLKcatPath);

% First check if the Docker image exist
[checks.image.status, checks.image.out] = system('docker images');
if checks.image.status == 0 && ~contains(checks.image.out,'dlkcat')
    % Check if the image file is already dowloaded
    if ~exist(fullfile(DLKcatPath,'dlkcat_docker.tar.gz'),'file')
        disp('=== Downloading DLKcat docker image, this may take more than 20 minutes, depending on your internet connection speed...')
        packageURL = 'https://github.com/SysBioChalmers/GECKO/raw/dlkcatPackage/dlkcat_docker.tar.gz';
        websave(fullfile(DLKcatPath,'dlkcat_docker.tar.gz'),packageURL,webOptions);
    end
    [checks.dockerload.status, checks.dockerload.out] = system('docker load --input dlkcat_docker.tar.gz');
    if deleteImage
        delete(fullfile(DLKcatPath,'dlkcat_docker.tar.gz'));
    end
    if checks.dockerload.status ~= 0
        error('Fail to load the image')
    end
elseif checks.image.status ~= 0
    error('Looks like Docker is installed, but make ensure that the docker desktop app is running.')
end

% Copy the input file for needed for the docker image to run DLKcat
copyfile(DLKcatFile,fullfile(DLKcatPath, 'data/DLKcat.tsv'));

disp('=== Running DLKcat prediction, this may take several minutes...')
% In the next line, pythonPath does not need to be specified, because it is
% already mentioned when building the virtualenv.
dlkcat.status = system('docker compose -f docker-compose.yaml up dlkcat-cont', '-echo');

if dlkcat.status == 0
    movefile('data/DLKcatOutput.tsv',DLKcatFile);
    cd(currPath);
else
    cd(currPath);
    error('DLKcat encountered an error.')
end

end
