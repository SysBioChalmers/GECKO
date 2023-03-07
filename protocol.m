%% This file accompanies the GECKO3 protocol (DOI TO BE ADDED)
% Where the commands shown in the protocol paper refer to generic
% functionality, these functions are used here to reconstruct, simulate and
% analyze ecYeastGEM, the ecModel of Saccharomyces cerevisiae.
%
% NOTE: this file might have more up-to-date information about the
% functions, and should always work with the most recent GECKO3 release.
% This file will have additional commands and analysis, that are not
% included as such in the protocol.

%% Preparation stage for ecModel reconstruction
% - Install GECKO % RAVEN
%   - GECKO can be installed via cloning or direct download of ZIP file.
%     See installation instructions #here (TO BE ADDED)
%   - Add the appropriate GECKO (sub)folders to MATLAB path:
GECKOInstaller.install %
%   - Install RAVEN by following the installation instructions:
%     https://github.com/SysBioChalmers/RAVEN/wiki/Installation
%   - The installation of Gurobi as LP solver is highly recommended
checkInstallation % Confirm that RAVEN is functional, should be 2.7.12 or later.

% - Prepare a folder in userData, or select another path whether all model-
%   specific files and scripts will be kept. This is already done for
%   ecYeastGEM.
% - Find a high-quality GEM of your species of interest. ecYeastGEM is
%   based on yeast-GEM https://github.com/SysBioChalmers/yeast-GEM/releases
% - Release v8.6.2 of yeast-GEM is also distributed with GECKO at 
%   userData/ecYeastGEM/model/yeast-GEM.xml
% - Modify the model adapter (at userData/ecYeastGEM/ecYeastGEMadapter.m)
%   to contain organism- and model-specific parameters.

%% STAGE 1: expansion from a starting metabolic model to an ecModel structure

% STEP 1 Set modelAdapter
geckoRoot = findGECKOroot;
ModelAdapterManager.setDefaultAdapterFromPath(fullfile(geckoRoot,'userData','ecYeastGEM')); 
% If you explicitly want to set the ModelAdapter as input parameter for
% many of the GECKO functions, you can get the ModelAdapter structure with
% the following command that should be run after the above command.
ModelAdapter = ModelAdapterManager.getDefaultAdapter();
% Make it easy to check what value the parameters are.
params = ModelAdapter.getParameters();

% STEP 2 Load conventional yeast-GEM
% If the location to the conventional GEM was already set in the modelAdapter,
% as obj.param.convGEM, then loadConventionalGEM can directly be used. If
% the model is stored somewhere else (and not specified in obj.param.convGEM),
% you can also use RAVEN's importModel(). In that case you will never use
% loadConventionalGEM and the obj.param.convGEM never has to be specified.
modelY = loadConventionalGEM();
% modelY = importModel(fullfile(modelRoot,'models','yeast-GEM.xml')); %Alternative

% STEP 3 Prepare ecModel
% We will make a full and a light GECKO model
[ecModelFull, noUniprot] = makeEcModel(modelY,false);
% Note that noUniprot is empty: for all genes a match could be find in the
% Uniprot dataset
[ecModelLight, noUniprot] = makeEcModel(modelY,true);
% Whil we initiate a light model here, for the remainder of this protocol.m
% file we will continue with the full version
ecModel = ecModelFull;
doc makeEcModel
% Read makeEcModel documentation to get a list of all it does: it prepare
% the new model.ec structure and prepares the S-matrix by splitting
% reversible reactions, isozymes etc.

% STEP 4 Annotate with complex data
% As Saccharomyces cerevisiae is available on Complex Portal
% (https://www.ebi.ac.uk/complexportal/), its data is included in
% ecYeastGEM.
complexInfo = getComplexData(); % No need to run, as data is already in adapter folder
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel, complexInfo);

% STEP 5 Store model in YAML format
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');
ecModel=loadEcModel('ecYeastGEMfull.yml');

%% STAGE 2: integration of kcat into the ecModel structure
% STEP 6 Decide which kcat source to use
% In the steps below, all options are shown. Not all are required, it is up
% to the user to decide which ones they want to use.

% STEP 7 Fuzzy matching with BRENDA
% Requires EC numbers, which are here first taken from the starting model,
% with the missing ones taken from Uniprot & KEGG databases.
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel);
% Do the actual fuzzy matching with BRENDA
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);
% Now we have a kcatList, which will be used to update ecModel in a later
% step.

% STEP 8 DLKcat prediction through machine learning
% Requires metabolite SMILES, which are gathered from PubChem
ecModel = findMetSmiles(ecModel);
% DLKcat runs in Python. An input file is written, which is then used by
% DLKcat, while the output file is read back into MATLAB.
writeDLKcatInput(ecModel);
% runDLKcat will attempt to install and run DLKcat, but this might not work
% for all systems. In that case, the user should install Docker Desktop
% (https://docs.docker.com/get-docker/) and use runDLKcatFromDocker()
% instead.
runDLKcat();
%runDLKcatFromDocker(); % Only uncomment and run if required.
kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 9 Combine kcat from BRENDA and DLKcat
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

% STEP 10 Take kcatList and populate edModel.ec.kcat
ecModel  = selectKcatValue(ecModel, kcatList_merged);

% STEP 11 Apply custom kcat values
% During the development of yeast-GEM ecModels (through GECKO 1 & 2),
% various kcat values have been manually curated, to result in a model that
% is able to validate a wide range of phenotypes. These curations are
% summarized under /data/customKcats.tsv, and applied here.
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);
% To modify the S-matrix, to actually implement the kcat/MW constraints,
% you run applyKcatConstraints.
ecModel  = applyKcatConstraints(ecModel);

% STEP 12 Get kcat values across isoenzymes
ecModel = getKcatAcrossIsoenzymes(ecModel);

% STEP 13 Get standard kcat
% Assign an enzyme cost to reactions without gene assocation (except
% exchange, transport and pseudoreactions)
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% STEP 14 Apply kcat constraints from ecModel.ec.kcats to ecModel.S
% This function can be run at any point to re-apply the kcat contraints on
% the model. It also considers the complex data 
ecModel = applyKcatConstraints(ecModel);

% STEP 15 Set upper bound of protein pool
ecModel = setProtPoolSize(ecModel);

%% STAGE 3: model tuning
% Test whether the model is able to reach maximum growth if glucose uptake
% is unlimited. First set glucose uptake unconstraint
ecModel = setParam(ecModel,'lb','r_1714',-1000);
% Run FBA
sol = solveLP(ecModel,1)
% It reaches growth rate 0.0936, while it should be able to reach 0.41.
% We can also look at the exchange fluxes, but it does not inform use too
% much at this point. Interesting to see that there is quite some ethanol
% fermentation going on.
printFluxes(ecModel, sol.x)

% STEP 16 Relax protein pool constraint
protPoolIdx = strcmp(ecModel.rxns, 'prot_pool_exchange');
ecModel.ub(protPoolIdx) = 1000;
% Important to perform parsimonious FBA by setting minFlux to 1:
sol = solveLP(ecModel,1);
ecModel.ub(protPoolIdx) = sol.x(protPoolIdx);

% STEP 17 Sensitivity tuning
% First reset the protein pool constraint to a more realistic value,
% reverting STEP 16.
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

% STEP 18 Bayesian tuning
% Alternative to STEP 17.
%ecModel = bayesianSensitivityTuning(ecModel);

%% STAGE 4 integration of proteomics data into the ecModel
% STEP 19 Load proteomics data and constrain ecModel 
protData = loadProtData(3); %Number of replicates
ecModel = fillProtConcs(ecModel,protData);
ecModel = constrainProtConcs(ecModel);

% STEP 20 Update protein pool
ecModel = updateProtPool(ecModel);

% STEP 21 Load flux data
fluxData = loadFluxData();
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose'); % Use first condition
sol = solveLP(ecModel); % To observe if growth was reached
disp(['Growth rate that is reached: ' num2str(abs(sol.f))])

% Growth rate of 0.1 is by far not reached, flexibilize protein
% concentrations
[ecModel, flexProt] = flexibilizeProtConcs(ecModel,0.1,10);

% It gets stuck at 0.0889. It seems like protein abundances are not
% preventing the model to reach 0.1. First look if the conventional GEM is
% able to reach 0.1 with the same constraints:

modelY = constrainFluxData(modelY,fluxData);
sol = solveLP(modelY)

% It also only reaches 0.0889! So the metabolic network would not be able
% to adhere to all measured constraints. Perhaps there is something
% incorrect with the measurements? For now, we will limit the growth rate
% to 0.08885 in our search:
[ecModel, flexProt] = flexibilizeProtConcs(ecModel,0.08885,10);

% Growth is reached! Let's make sure we store this functional model
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');

%% STAGE 5: simulation and analysis
% If starting from here, load some basic assets
ecModel = loadEcModel('ecYeastGEMfull.yml',ModelAdapter);
modelY = loadConventionalGEM();
fluxData = loadFluxData;

% STEP 23 Example of various useful RAVEN functions
% % Set the upper bound of reaction r_0001 to 10.
% ecModel = setParam(ecModel,'ub','r_0001',10);
% % Set the lower bound of reaction r_0001 to 0.
% ecModel = setParam(ecModel,'lb','r_0001',0);
% % Set the objective function to maximize reaction biomassRxn
% ecModel = setParam(ecModel,'obj','r_4041',1);
% % Set the objective function to minimize protein usage
% ecModel = setParam(ecModel,'obj','prot_pool_exchange',-1);
% % Perform flux balance analysis (FBA)
% sol = solveLP(ecModel);
% % Perform parsimonious FBA (minimum total flux)
% sol = solveLP(ecModel,1);
% % Inspect exchange fluxes from FBA solution
% printFluxes(ecModel,sol.x)
% % Inspect all (non-zero) fluxes from FBA solution
% printFluxes(ecModel,sol.x,false)
% % Export to Excel file (will not contain potential model.ec content)
% exportToExcelFormat(ecModel,'filename.xlsx');

% STEP 24 Selecting objective functions
ecModel = setParam(ecModel,'obj',params.bioRxn,1);
sol = solveLP(ecModel)
disp(['Growth rate reached: ' num2str(abs(sol.f))])
% Set growth lower bound to 99% of the previous value
ecModel = setParam(ecModel,'lb',params.bioRxn,0.99*abs(sol.f));
ecModel = setParam(ecModel,'obj','prot_pool_exchange',-1);
sol = solveLP(ecModel)
disp(['Minimum protein pool usage: ' num2str(abs(sol.f)) ' ug/gDCW'])

% STEP 25 Compare fluxes from ecModel and starting model
% Constrain with the same conditions to model and ecModel. We now fix the
% observed growth as lower bound ('min' in constrainFluxData) and allow 5%
% variability around the other measured fluxes.
% We know that growth can only reach 0.088, so use this instead of 0.1.
fluxData.grRate(1) = 0.088;
ecModel = constrainFluxData(ecModel,fluxData,1,'min',5);
% Minimize protein pool usage.
ecModel = setParam(ecModel,'obj','prot_pool_exchange',-1);
solEC = solveLP(ecModel,1)

% Apply (almost) the same to non-ecModel. Same constraints on fluxes, but
% objective function remains growth (cannot minimize protein usage).
modelY = constrainFluxData(modelY,fluxData,1,'min',5);
sol = solveLP(modelY,1)
% Note that with the 5% of variability, the model can reach higher growth
% rate than 0.0885, as it is allowed to take up a little bit more glucose
% this time.

% Map the ecModel fluxes back to the non-ecModel
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, modelY, solEC.x);

% Print the fluxes next to each other, showing only exchange reactions
printFluxes(modelY,[sol.x mappedFlux])
% Print the fluxes next to each other, showing all reactions
printFluxes(modelY,[sol.x mappedFlux],false)
% Look at ratio of ecModel / non-ecModel
ratioFlux = mappedFlux ./ sol.x;
ratioFlux(isnan(ratioFlux)) = 0; % Divisions by zero give NaN, reset to zero.
printFluxes(modelY,ratioFlux,false)

% STEP 26 Inspect enzyme usage
% Show the result from the earlier simulation, without mapping to
% non-ecModel.
usageData = enzymeUsage(ecModel,solEC.x);
usageReport = reportEnzymeUsage(ecModel,usageData,0.90);

% STEP 27 Perform (ec)FVA
% Perform FVA
[minFluxEc, maxFluxEc] = ecFVA(ecModel, modelY);
[minFluxY, maxFluxY] = ecFVA(modelY, modelY);

% Write results to output file
output = [modelY.rxns, modelY.rxnNames, num2cell([minFluxY, maxFluxY, minFlux, maxFlux])]';
fID = fopen(fullfile(params.path,'output','ecFVA.tsv'),'w');
fprintf(fID,'%s %s %s %s %s %s\n','rxnIDs', 'rxnNames', 'minFlux', ...
            'maxFlux', 'ec-minFlux', 'ec-maxFlux');
fprintf(fID,'%s %s %g %g %g %g\n',output{:});
fclose(fID);

% Look at flux ranges to indicate reaction-level variability
fluxRange = maxFlux - minFlux;
fluxRangeY = maxFluxY - minFluxY;

% Plot variability distributions of both models in 1 plot
hold on
cdfplot(fluxRange)
cdfplot(fluxRangeY)
set(gca, 'XScale', 'log', 'Xlim', [1e-8 1e4])
set(findall(gca, 'Type', 'Line'), 'LineWidth', 2)
legend(['ecYeastGEM (mean: ' num2str(mean(fluxRange)) ')'],...
        ['yeast-GEM (mean: ' num2str(mean(fluxRangeY)) ')'],...
        'Location','northwest')
title('Flux variability (cumulative distribution)');
xlabel('Variability range [mmol/gDw h]');
ylabel('Cumulative distribution');
saveas(gca, fullfile(params.path,'output','ecFVA.pdf'))
