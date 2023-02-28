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
% The pipeline can look different dependent on user preferences. Getting an
% all-DLKcat model would look something like:
% makeEcModel -> findMetSmiles -> writeDLKcatInput -> readDLKcatOutput ->
% selectKcatValue -> applyKcatConstraints -> ...
% while a fuzzy matching to BRENDA approach looks like:
% makeEcModel -> getECfromDatabase -> fuzzyMatching -> applyKcatConstraints
% -> ...

%% Initiate model reconstruction

% Set the ModelAdapter correctly. This loads the ModelAdapter file that is
% in userData/ecYeastGEM/.
ModelAdapterManager.setDefaultAdapterFromPath(fullfile(findGECKOroot,'userData','ecYeastGEM'), 'true'); 
% If you explicitly want to set the ModelAdapter as input parameter for
% many of the GECKO functions, you can get the ModelAdapter structure with
% the following command that should be run after the first.
ModelAdapter = ModelAdapterManager.getDefaultAdapter();

% If the location to the conventional GEM was already set in the modelAdapter,
% as obj.param.convGEM, then loadConventionalGEM can directly be used. If
% the model is stored somewhere else (and not specified in obj.param.convGEM),
% you can also use RAVEN's importModel(). In that case you will never use
% loadConventionalGEM and the obj.param.convGEM never has to be specified.
modelY = loadConventionalGEM();
% modelY = importModel(fullfile(modelRoot,'models','yeast-GEM.xml')); %Alternative

% Prepare ec-model
[ecModel, noUniprot] = makeEcModel(modelY,false,ModelAdapter);
% Read makeEcModel documentation to get a list of all it does: it prepare
% the new model.ec structure and prepares the S-matrix by splitting
% reversible reactions, isozymes etc.

% Can also make a geckoLight model with makeEcModel(..,..,true);

% Store model in YAML format. Requires RAVEN 2.7.12. saveEcModel and
% loadEcModel take the obj.param.path, and assume that the model will be in
% the models subfolder, and named ecModel.yml. Alternatives names can be
% provided. RAVEN's readYAMLmodel and writeYAMLmodel can also be used, but
% then require specifying the whole path.
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEM');
ecModel=loadEcModel('ecYeastGEM.yml');

%% Gather complex data
% For species with data in ComplexPortal, you can gather that information
% with getComplexData (no need to repeat for yeast-GEM, already stored in
% data subfolder)
% complexInfo = getComplexData();
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);


%% Gather kcat values
% Different approaches are possible: (1) DLKcat; (2) fuzzy matching; (3)
% manual kcatList. For (1) and (2), ecRxns parameter can indicate that you
% do not want to get kcats for all reactions. You might want to first run
% fuzzy matching, only keep the ones without any EC wildcards, and then run
% DLKcat for the remaining reactions.

% (1) DLKcat
% Requires metabolite SMILES:
ecModel = findMetSmiles(ecModel);

% A DLKcat.tsv file is written, which is later populated by DLKcat with
% kcat values. If the file already exists, it will not be overwritten, to
% avoid losing existing kcat values (unless 'overwrite' was set as 'true'
% when running writeDLKcatInput.
% runDLKcat will attempt to download, install and run DLKcat, but this
% might not work for all systems. In that case, the user will be directed
% to manually download, install and DLKcat via the GECKO-provided DLKcat package

writeDLKcatInput(ecModel);
runDLKcat();
kcatList_DLKcat = readDLKcatOutput(ecModel);

% (2) fuzzy matching
% Requires EC-codes. Can be either gathered from the model.eccodes field,
% or from the Uniprot (and KEGG) data
ecModel         = getECfromGEM(ecModel);
ecModel         = getECfromDatabase(ecModel);
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% (3) combine fuzzy matching and DLkcat
% Assumes that you've run both step (1) and step (2)
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);
% In kcat list, multiple kcat values might be provided for each reaction,
% for instance because of multiple substrates (as coming from DLKcat).
% selectKcatValue selects one kcat value per reaction and populates
% model.ec.kcat.
ecModel_merged  = selectKcatValue(ecModel, kcatList_merged);
% To modify the S-matrix, to actually implement the kcat/MW constraints,
% you run applyKcatConstraints.
ecModel_merged  = applyKcatConstraints(ecModel_merged);

% Additional kcats. In previous work, kcat values of various reactions were
% manual curated. Here we will apply those custom kcats.
ecModel_merged = applyCustomKcats(ecModel_merged); % Populate model.ec.kcats
ecModel_merged = applyKcatConstraints(ecModel_merged);

% Sometimes isoenzymes are not assigned kcat values. This is resolved here
ecModel_merged = getKcatAcrossIsoenzymes(ecModel_merged);

% Reactions without gene association have no protein cost. To counter this,
% a standard kcat and MW are assigned to those reactions (except for e.g.
% exchange and pseudoreactions).
ecModel_merged = getStandardKcat(ecModel_merged);
% Re-apply the updated model.ec.kcats
ecModel_merged = applyKcatConstraints(ecModel_merged);

%% Do some first simulations
% Set glucose unlimited
ecModel_merged = setParam(ecModel_merged,'lb','r_1714',-1000);
% Set total protein constraint
ecModel_merged = setProtPoolSize(ecModel_merged);

sol = solveLP(ecModel_merged,1)
printFluxes(ecModel_merged, sol.x)
% Growth rate is not high enough (0.0989 instead of 0.41).
%% Tune kcat values to reach max growth rate
% Let's increase the kcat values of reactions with the most-used enzymes by
% 2-fold in each iteration. Most-used enzymes (% of protein pool) are most-
% limiting growth. A reaction might have its kcat increased multiple times.
% TunedKcats gives an overview of what kcat values were changed
[ecModelTuned, tunedKcats] = sensitivityTuning(ecModel_merged);

% Tune sigma-factor (not sure why this is done, as it is set to 0.5 default
% and it just reaches the same here again?
[ecModelTuned, optSigma] = sigmaFitter(ecModelTuned);

%% Contrain with proteomics data
% Load proteomics
protData = loadProtData(3); %Number of replicates
ecModelProt = fillProtConcs(ecModelTuned,protData);
ecModelProt = updateProtPool(ecModelProt);
ecModelProt = constrainProtConcs(ecModelProt);

% Load matching flux data
fluxData = loadFluxData();
ecModelProtFlux = constrainFluxData(ecModelProt,fluxData);
sol = solveLP(ecModelProtFlux)
% Growth rate of 0.1 is by far not reached, flexibilize protein
% concentrations
[ecModelFlex, flexProt] = flexibilizeProtConcs(ecModelProtFlux,0.1,10);

% It gets stuck at 0.0889. It seems like protein abundances are not
% preventing the model to reach 0.1. First look if the conventional GEM is
% able to reach 0.1 with the same constraints:
modelY = constrainFluxData(modelY,fluxData);
sol = solveLP(modelY)
% It also only reaches 0.0889! So the metabolic network would not be able
% to adhere to all measured constraints. Perhaps there is something
% incorrect with the measurements? For now, we will limit the growth rate
% to 0.08885 in our search:
[ecModelFlex, flexProt] = flexibilizeProtConcs(ecModelProtFlux,0.08885,10);

% Growth is reached!
%ecModelProtFlux = bayesianSensitivityTuning(ecModelProtFlux);

%% Perform simulations
% Constrain with the same conditions to model and ecModel
% We know that growth can only reach 0.088
fluxData.grRate(1) = 0.088;
ecModelProtFlux = constrainFluxData(ecModelFlex,fluxData,1,'min',5);
solEC = solveLP(ecModelProtFlux,1)
modelY = constrainFluxData(modelY,fluxData,1,'min',5);
sol = solveLP(modelY,1)
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModelProtFlux, modelY, solEC.x);

printFluxes(modelY,[sol.x mappedFlux])

% Perform FVA
[minFlux, maxFlux] = ecFVA(ecModelProtFlux, modelY);

