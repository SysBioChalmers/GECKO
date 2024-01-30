function runDLKcat(modelAdapter)
% runDLKcat
%   Runs DLKcat to predict kcat values from a Docker image. Once DLKcat is succesfully
%   run, the DLKcatFile will be overwritten with the DLKcat
%   output in the model-specific 'data' sub-folder taken from modelAdapter
%   (e.g. GECKO/tutorials/tutorial_yeast-GEM/data/DLKcat.tsv)
%
% Input
%   modelAdapter    a loaded model adapter. (Optional, will otherwise use
%                   the default model adapter)
%
%   NOTE: 1. Requires Docker to be installed, and Docker Desktop running. Visit "https://www.docker.com"
%         2. Runtime will depend on whether the image is to be downloaded or not.

if nargin < 1 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

params = modelAdapter.params;

%% Check and install requirements
% On macOS, Docker might not be properly loaded if MATLAB is started via
% launcher and not terminal.
if ismac
    setenv('PATH', strcat('/usr/local/bin', ':', getenv("PATH")));
end

% Check if Docker is installed
[checks.docker.status, checks.docker.out] = system('docker --version');
if checks.docker.status ~= 0
    error('Cannot find Docker, make sure it is installed. If it is, it might be required to start Matlab from the command-line instead of the launcher in order for Docker to be detected and used.')
end

disp('Running DLKcat prediction, this may take many minutes, especially the first time.')
status = system(['docker run --rm -v "' fullfile(params.path,'/data') '":/data ghcr.io/sysbiochalmers/dlkcat-gecko:0.1 /bin/bash -c "python DLKcat.py /data/DLKcat.tsv /data/DLKcatOutput.tsv"']);

if status == 0 && exist(fullfile(params.path,'data/DLKcatOutput.tsv'))
    delete(fullfile(params.path,'/data/DLKcat.tsv'));
    movefile(fullfile(params.path,'/data/DLKcatOutput.tsv'), fullfile(params.path,'/data/DLKcat.tsv'));
    disp('DKLcat prediction completed.');
else    
    error('DLKcat encountered an error or it did not create any output file.')
end

