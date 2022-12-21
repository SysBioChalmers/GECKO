function [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,parameters,name,modelVer)
% enhanceGEM
%
%   Main function for running the GECKO pipeline. It returns an ecModel and
%   its constrained version with an upper limit in the total protein pool
%   usage pseudoreaction (ecModel_batch) with calibrated Kcat parameters that
%   allow it to grow at a specified experimental growth rate.
%
%   model       a GEM MATLAB structure compatible with the COBRA or RAVEN
%               toolbox.
%   toolbox     string with the name of the prefered toolbox for model SBML
%               export (COBRA or RAVEN)
%   name        Desired name for the ecModel (opt, default '')
%   modelVer    modelVer of the original GEM (opt, default '')
%
%
%   ecModel        an ecModel MATLAB structure suitable for incorporation of   
%                  proteomics data as individual enzyme usage constraints.
%   ecModel_batch  an ecModel MATLAB structure with a global constraint on 
%                  the total protein pool usage pseudoreaction,
%                  proportional to the measured total protein content (Ptot)
%
%   Usage: [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,modelVer)
%
%   Ivan Domenzain. Last edited: 2020-10-05
%

current = pwd;

% Save geckomat path as a parameter to return after use organism specific
% scripts
parameters.geckomat = current;

if nargin < 3
    error('Not enough input arguments, parameters are required');
end

if nargin < 4
    name    = 'ecModel';
end

if nargin < 5
    if isfield(model,'id')
        modelVer = model.id;
    else
        modelVer = '1';
    end
end

% Validtae if there is a defined path where organism/model specific scripts are
if ~isfield(parameters, 'userDataPath')
    [~,values] = fileattrib('../../userData/');
    % Add to parameters the full path
    parameters.userDataPath = [values.Name '/' name '/'];
else
    [~,values] = fileattrib(parameters.userDataPath);
    %Update to a full path
    parameters.userDataPath = [values.Name '/' name '/'];
end

% Validtae if there is a defined path where save the output
if ~isfield(parameters, 'outputPath')
    [~,values] = fileattrib('../../userData/');
    % Add to parameters the full path
    parameters.outputPath = [values.Name '/' name '/output/'];
else
    [~,values] = fileattrib(parameters.outputPath);
    %Update to a full path
    parameters.outputPath = [values.Name '/' name '/'];
end

%create specific subfolder for ecModel files if not present already
if ~isfolder(parameters.outputPath)
    mkdir(parameters.outputPath)
end

%Convert model to RAVEN for easier visualization later on:
format short e
if isfield(model,'rules')
    initCobraToolbox
    model = ravenCobraWrapper(model);
end

fprintf('\n***************************************************************')
fprintf('\n   GECKO: Adding enzyme constraints to a genome-scale model')
fprintf('\n***************************************************************\n\n')

%Remove blocked rxns + correct model.rev:
cd([current '/change_model'])
[model,~,~] = preprocessModel(model,name,modelVer);

fprintf('\n==================')
fprintf('\nGenerating ecModel:')
fprintf('\n==================\n')

%Retrieve kcats & MWs for each rxn in model:
cd([current '/get_enzyme_data'])
model_data = getEnzymeCodes(model,parameters.userDataPath);
kcats      = matchKcats(model_data,parameters.org_name);

%Integrate enzymes in the model:
cd([current '/change_model'])
ecModel                 = readKcatData(model_data,kcats,parameters);
cd([parameters.userDataPath 'scripts/'])
[ecModel,modifications] = manualModifications(ecModel,parameters);

%Constrain model to batch conditions:
fprintf('\n==============================================================')
fprintf('\nGenerating ecModel with shared pool assumption (ecModel_batch):')
fprintf('\n==============================================================\n')
cd([current '/limit_proteins'])
[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,modifications,name,parameters);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

%Save output models:
fprintf('\n=============')
fprintf('\nSaving models:')
fprintf('\n=============\n')
cd(current)
ecModel = saveECmodel(ecModel,toolbox,name,modelVer,parameters.outputPath);
ecModel_batch = saveECmodel(ecModel_batch,toolbox,[name '_batch'],modelVer,parameters.outputPath);

end