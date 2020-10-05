function [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,version)
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
%   version     version of the original GEM (opt, default '')
%
%
%   ecModel        an ecModel MATLAB structure suitable for incorporation of   
%                  proteomics data as individual enzyme usage constraints.
%   ecModel_batch  an ecModel MATLAB structure with a global constraint on 
%                  the total protein pool usage pseudoreaction,
%                  proportional to the measured total protein content (Ptot)
%
%   Usage: [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,version)
%
%   Ivan Domenzain. Last edited: 2019-07-13
%

if nargin < 3
    name    = '';
end
if nargin < 4
    version = '';
end

%Convert model to RAVEN for easier visualization later on:
format short e
if isfield(model,'rules')
    initCobraToolbox
    model = ravenCobraWrapper(model);
end
%Get model-specific parameters
parameters = getModelParameters;
%Remove blocked rxns + correct model.rev:
cd change_model
[model,name,version] = preprocessModel(model,name,version);

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,parameters.org_name);

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);

%create specific subfolder for ecModel files if not present already
if ~isdir(['../models/' name])
    mkdir(['../models/' name])
end
%Constrain model to batch conditions:
cd ../limit_proteins
[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

%Save output models:
cd ../../models
ecModel = saveECmodel(ecModel,toolbox,name,version);
ecModel_batch = saveECmodel(ecModel_batch,toolbox,[name '_batch'],version);
cd ../geckomat

end