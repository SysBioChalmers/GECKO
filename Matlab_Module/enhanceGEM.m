%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name,version)
%
% Benjamin J. Sanchez & Ivan Domenzain. Last edited: 2018-03-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name,version)

%Provide your organism scientific name
org_name = 'saccharomyces cerevisiae';
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
kcats      = matchKcats(model_data,org_name);
save(['../../Models/' name '/data/' name '_enzData.mat'],'model_data','kcats','version')

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);

%Constrain model to batch conditions:
sigma  = 0.5;      %Optimized for glucose
Ptot   = 0.5;      %Assumed constant
gR_exp = 0.41;     %[g/gDw h] Max batch gRate on minimal glucose media
cd ../limit_proteins
[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,sigma,Ptot,gR_exp,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

%Save output models:
cd ../../Models
ecModel.description       = [name '_' version];
ecModel_batch.description = [name '_batch_' version];
save([name '/' name '.mat'],'ecModel')
save([name '/' name '_batch.mat'],'ecModel_batch')
saveECmodelSBML(ecModel,name,false);
saveECmodelSBML(ecModel_batch,name,true);
cd ../Matlab_Module

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%