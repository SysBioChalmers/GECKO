%% This script sketches the blueprint for the GECKO3 protocol.
% At a later stage, it can likely be converted into a script that is
% referenced in the protocol paper. Right now, this script might contain
% comments and clarifications that should rather be included in the
% protocol than in here, while other comments might be best kept in this
% script. At the moment it is just a sketchbook.

%% Prepare for ec model reconstruction
% - Prepare a folder in userData, or select another path whether all model-
%   specific files and scripts will be kept
% - Find a high-quality GEM of your species of interest
%   - For reference, might be best to be stored together with the
%     model-specific files, e.g. userData/ecYeastGEM/model/yeast-GEM.xml
% - Install RAVEN & GECKO
%   - RAVEN can be installed in many different ways, refer to GitHub Wiki.
checkInstallation % To confirm that RAVEN is functioning

%   - GECKO can be installed via cloning or direct download of ZIP file.
%     Should GECKO also be distributed as MATLAB Add-On, just like RAVEN?
%   - GECKO should be added to the MATLAB path.
GECKOInstaller.install % Adds the appropriate folders to MATLAB

% - On Uniprot:
%   - Search a proteome for your species and note the proteome ID
%   - Check which gene identifiers are used in the model, and see where
%     this is represented in Uniprot. E.g. for yeast, the genes are in
%     "Gene Names (Ordered locus)", which via the API can be gather via
%     "gene_oln" according to https://www.uniprot.org/help/return_fields
% - On KEGG:
%   - Find the organism code of your species
% - Modify the ModelAdapter file in the appropriate userData subfolder 
%   with some of the above information and more.

%% Summary
% Pipeline can look different dependent on user preferences. Getting an
% all-DLKcat model would look something like:
% makeEcModel -> findMetSmiles -> writeDLKcatInput -> readDLKcatOutput ->
% selectKcatValue -> applyKcatConstraints -> ...
% while a fuzzy matching to BRENDA approach looks like:
% makeEcModel -> getECfromDatabase -> fuzzyMatching -> applyKcatConstraints
% -> ...

%% Initiate model reconstruction

% Set the ModelAdapter correctly. This loads the ModelAdapter file that is
% in userData/ecYeastGEM/. First define modelRoot as userData/ecYeastGEM
modelRoot = fullfile(findGECKOroot,'userData','ecYeastGEM');
ModelAdapterManager.setDefaultAdapterFromPath(fullfile(modelRoot), 'true'); 
ModelAdapter = ModelAdapterManager.getDefaultAdapter();

% Load the model with RAVEN's importModel
% For yeast-GEM, the model is already in userData/ecYeastGEM/models/
modelY = importModel(fullfile(modelRoot,'models','yeast-GEM.xml'));

% Prepare ec-model
ecModel = makeEcModel(modelY);
% Read makeEcModel documentation to get a list of all it does: it prepare
% the new model.ec structure and prepares the S-matrix by splitting
% reversible reactions, isozymes etc.

% Can also make a geckoLight model with makeEcModel(..,..,true);

% Store model in YAML format. Might be good to write a small wrapper
% function, so you do not need to provide the path from ModelAdapter.
% writeYAMLmodel is a RAVEN 2.7.10 function
writeYAMLmodel(ecModel,fullfile(ModelAdapter.params.path,'models','ecYeastGEM'));
ecModel=readYAMLmodel(fullfile(ModelAdapter.params.path,'models','ecYeastGEM.yml'));

%% Gather complex data
% For species with data in ComplexPortal, you can gather that information
% with getComplexData (no need to repeat for yeast-GEM, already stored in
% data subfolder)
% complexInfo = getComplexData();
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

% Describe manual curation. Albert has a function to also apply text file
% with curated complex data.

%% Gather kcat values
% Different approaches are possible: (1) DLKcat; (2) fuzzy matching; (3)
% manual kcatList. For (1) and (2), ecRxns parameter can indicate that you
% do not want to get kcats for all reactions. You might want to first run
% fuzzy matching, only keep the ones without any EC wildcards, and then run
% DLKcat for the remaining reactions.

% (1) DLKcat
% Requires metabolite SMILES:
ecModel = findMetSmiles(ecModel);

% Currently, a DLKcatInput.tsv file is written that can be used by DLKcat,
% and the DLKcatOutput.tsv file can be loaded into MATLAB again. runDLKcat
% will attempt to download, install and run DLKcat, but this might not work
% for all systems. In that case, the user will be directed to manually
% download, install and DLKcat via the GECKO-provided DLKcat package

writeDLKcatInput(ecModel);
runDLKcat();
kcatList_DLKcat = readDLKcatOutput(ecModel);
ecModel         = selectKcatValue(ecModel,kcatList_DLKcat);
ecModel         = applyKcatConstraints(ecModel);

% (2) fuzzy matching
% Requires EC-codes. Can be either gathered from the model.eccodes field,
% or from the Uniprot (and KEGG) data
%ecModel         = getECfromGEM(ecModel);
ecModel         = getECfromDatabase(ecModel);
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);
ecModel_fuzzy   = selectKcatValue(ecModel, kcatList_fuzzy);
ecModel_fuzzy   = applyKcatConstraints(ecModel_fuzzy);

% (3) combine fuzzy matching and DLkcat
% Assumes that you've run both step (1) and step (2)
kcatList_merged = mergeDlkcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);
ecModel_merged  = selectKcatValue(ecModel, kcatList_merged);
ecModel_merged  = applyKcatConstraints(ecModel_merged);

writeYAMLmodel(ecModel_merged,fullfile(ModelAdapter.params.path,'models','ecYeastGEM'));

% Solve without proteome constraint
sol = solveLP(ecModel_merged,1)
printFluxes(ecModel_merged, sol.x)
% RAVEN currently gives warning about _REV in reaction identifiers when ,1
% is used in solveLP, this will be removed in next RAVEN release

%% Tune kcat values to reach max growth rate
% Protein = 0.5; enzyme = 0.5; saturation = 0.5; = 0.125
ecModel_merged = setProtPoolSize(ecModel_merged);
% Unlimited glucose uptake
ecModel_merged = setParam(ecModel_merged,'lb','r_1714',-1000);
sol = solveLP(ecModel_merged)
% Growth rate is not high enough (0.01 instead of 0.41). Let increase the
% kcat values of reactions with the most-used enzymes by 2-fold in each
% iteration. Most-used enzymes (% of protein pool) are most-limiting growth.
% A reaction might have its kcat increased multiple times. TunedKcats gives
% an overview of what kcat values were changed
[ecModelTuned, tunedKcats] = sensitivityTuning(ecModel_merged,[],[],2);

% Tune sigma-factor (not sure why this is done, as it is set to 0.5 default
% and it just reaches the same here again?
[ecModelTuned, optSigma] = sigmaFitter(ecModelTuned);

%% Set realistic conditions
% Model-specific reactions/scripts can be defined for this, but this is so
% specific that it does not make much sense to make a generic function like
% changeMedia_yeast



%% Contrain with proteomics data
% Load proteomics
ecModelProt = readProteomics(ecModelTuned);
ecModelProt = constrainProtConcs(ecModelProt);
sol=solveLP(ecModelProt)