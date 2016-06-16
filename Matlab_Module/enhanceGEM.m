%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [enzymeModel,model_data,kcats] = enhanceGEM(model,toolbox)
% 
%
% Benjamín J. Sánchez. Last edited: 2015-09-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [enzymeModel,model_data,kcats] = enhanceGEM(model,toolbox)

format short e
if strcmp(toolbox,'COBRA')
    initCobraToolbox
end
cd get_enzyme_data

%Update yeast 7 model with all recommended changes:
model = modelCorrections(model);

%Add some RAVEN fields for easier visualization later on:
model = standardizeModel(model,toolbox);

%Retrieve kcats & MWs for each rxn in model:
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data);

%Integrate enzymes in the model:
cd ../change_model
enzymeModel = readKcatData(model_data,kcats);
enzymeModel = manualModifications(enzymeModel);
cd ..

%Finally, save output model:
save([enzymeModel.description(1:end-4) '_max_current.mat'],'enzymeModel','model_data','kcats')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%