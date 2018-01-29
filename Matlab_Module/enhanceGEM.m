%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name)
%
% Benjam?n J. S?nchez. Last edited: 2017-04-12
%Ivan Domenzain.       Last edited: 2018-01-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,org_name)
Protdatabase  = 'ProtDatabase.mat';
%Provide your organism scientific name
org_name      = 'saccharomyces cerevisiae';
org_code      = 'sce';
format short e

if strcmp(toolbox,'COBRA')
   initCobraToolbox
end

%Update yeast 7 model with all recommended changes:
cd get_enzyme_data
model = modelCorrections(model);

%Add some RAVEN fields for easier visualization later on:
model = standardizeModel(model,toolbox);

%Retrieve kcats & MWs for each rxn in model:
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data, org_name);
cd ../../Models
save([org_code '_enzData.mat'],'model_data','kcats')
%Integrate enzymes in the model:
cd ../Matlab_Module/change_model
ecModel = readKcatData(model_data,kcats);
ecModel = manualModifications(ecModel);

%Constrain model to batch conditions:
cd ../limit_proteins
sigma         = 0.51;      %Optimized for glucose
Ptot          = 0.5;       %Assumed constant
ecModel_batch = constrainEnzymes(ecModel,Ptot,sigma);

Save output models:
 cd ../../models
 save([name '.mat'],'ecModel','model_data','kcats')
 save([name '_batch.mat'],'ecModel_batch')
saveECmodelSBML(ecModel,name);
saveECmodelSBML(ecModel_batch,[name '_batch']);
cd ../Matlab_Module

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%