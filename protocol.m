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
% - Install GECKO & RAVEN
%   - GECKO can be installed via cloning or direct download of ZIP file.
%     See installation instructions in the README.md:
%     https://github.com/SysBioChalmers/GECKO/tree/main#installation
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
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);
% complexInfo can be given as second input, but not needed here, as it will
% read the file that was written by getComplexData

% STEP 5 Store model in YAML format
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');
ecModel=loadEcModel('ecYeastGEMfull.yml');

%% STAGE 2: integration of kcat into the ecModel structure
% Decide which kcat source to use. In the steps below, all options are
% shown. Not all are required, it is up to the user to decide which ones
% they want to use.

% STEP 6 Fuzzy matching with BRENDA
% Requires EC numbers, which are here first taken from the starting model,
% with the missing ones taken from Uniprot & KEGG databases.
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel,noEC);
% Do the actual fuzzy matching with BRENDA
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);
% Now we have a kcatList, which will be used to update ecModel in a later
% step.

% STEP 7 DLKcat prediction through machine learning
% Requires metabolite SMILES, which are gathered from PubChem
ecModel = findMetSmiles(ecModel);
% DLKcat runs in Python. An input file is written, which is then used by
% DLKcat, while the output file is read back into MATLAB.
writeDLKcatInput(ecModel);
% runDLKcat will run the DLKcat algorithm via a Docker image
runDLKcat();
kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 8 Combine kcat from BRENDA and DLKcat
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

% STEP 9 Take kcatList and populate edModel.ec.kcat
ecModel  = selectKcatValue(ecModel, kcatList_merged);

% STEP 10 Apply custom kcat values
% During the development of yeast-GEM ecModels (through GECKO 1 & 2),
% various kcat values have been manually curated, to result in a model that
% is able to validate a wide range of phenotypes. These curations are
% summarized under /data/customKcats.tsv, and applied here.
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);
% To modify the S-matrix, to actually implement the kcat/MW constraints,
% applyKcatConstraints should be run. This takes the values from
% ecModel.ec.kcat, considers any complex/subunit data that is tracked in
% ecModel.ec.rxnEnzMat, together with the MW in ecModel.ec.mw, and uses
% this to modify the enzyme usage coefficients directly in ecModel.S. Any
% time a change is made to the .kcat, .rxnEnzMat or .mw fields, the
% applyKcatConstraints function should be run again to reapply the new
% constraints onto the metabolic model.
ecModel  = applyKcatConstraints(ecModel);

% STEP 11 Get kcat values across isoenzymes
ecModel = getKcatAcrossIsoenzymes(ecModel);

% STEP 12 Get standard kcat
% Assign an enzyme cost to reactions without gene assocation (except
% exchange, transport and pseudoreactions). These reactions are identified
% as those with empty entry in ecModel.grRules
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
% As the above functions have modified ecModel.ec.kcat,
% applyKcatConstraints is rerun as explained in step 11.
ecModel = applyKcatConstraints(ecModel);

% STEP 14 Set upper bound of protein pool
% The protein pool exchange is constrained by the total protein content
% (Ptot), multiplied by the f-factor (ratio of enzymes/proteins) and the
% sigma-factor (how saturated enzymes are on average: how close to their
% Vmax to they function based on e.g. metabolite concentrations). In 
% modelAdapter Ptot, f- and sigma-factors can all be specified (as rough
% estimates, 0.5 for each of the three parameters is reasonable).
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;
% But these values can also be defined separately. The f-factor can be 
% calculated from quantitative proteomics data, for instance with data that
% is available via PAXdb (https://pax-db.org/). calculateFfactor can be used to estimate the f-factor.
%f = calculateFfactor(ecModel); % Optional

ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

% Note that at a later stage (after stage 3), the sigma factor be further
% adjusted with sigmaFitter, to get a model that is able to reach a
% particular maximum growth rate. This will not be done here, as we first
% need to tune the kcat values in Stage 3.

%% STAGE 3: model tuning
% Test whether the model is able to reach maximum growth if glucose uptake
% is unlimited. First set glucose uptake unconstraint
ecModel = setParam(ecModel,'lb','r_1714',-1000);
% Run FBA
sol = solveLP(ecModel,1)
% It reaches growth rate 0.0877, while it should be able to reach 0.41 (the
% maximum growth rate of S. cerevisiae, that is entered in the model
% adapter as obj.params.gR_exp. We can also look at the exchange fluxes,
% but it does not inform use too much at this point. Interesting to see
% that there is quite some ethanol fermentation going on.
printFluxes(ecModel, sol.x)

% STEP 15 Relax protein pool constraint
% As a simplistic way to ensure the model to reach the growth rate, the
% upper bound of the protein pool exchange reaction can be increased to
% whatever is required. This works, but STEP 16 is preferred.
protPoolIdx = strcmp(ecModel.rxns, 'prot_pool_exchange');
ecModel.lb(protPoolIdx) = -1000;
% Important to perform parsimonious FBA by setting minFlux to 1:
sol = solveLP(ecModel,1)
ecModel.lb(protPoolIdx) = sol.x(protPoolIdx);

% STEP 16 Sensitivity tuning
% First reset the protein pool constraint to a more realistic value,
% reverting STEP 16.
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

%% STAGE 4 integration of proteomics data into the ecModel
% STEP 17 Load proteomics data and constrain ecModel 
protData = loadProtData(3); %Number of replicates
ecModel = fillProtConcs(ecModel,protData);
ecModel = constrainProtConcs(ecModel);

% STEP 18 Update protein pool
ecModel = updateProtPool(ecModel);

% STEP 19 Load flux data
fluxData = loadFluxData();
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose'); % Use first condition
sol = solveLP(ecModel); % To observe if growth was reached
disp(['Growth rate that is reached: ' num2str(abs(sol.f))])
% Growth rate of 0.1 is by far not reached, flexibilize protein
% concentrations

% STEP 20 Protein concentrations are flexibilized (increased), until the
% intended growth rate is reached.
[ecModel, flexProt] = flexibilizeProtConcs(ecModel,0.1,10);

% Neither individual protein levels nor total protein pool are limiting
% growth. Test whether the starting model is able to reach 0.1.
modelY = constrainFluxData(modelY,fluxData);
sol = solveLP(modelY)

% It also only reaches 0.0889! So the metabolic network would not be able
% to adhere to all measured constraints. Perhaps there is something
% incorrect with the measurements? Regardless, the ecModel is now able to
% reach about 0.0889, which will be fine for now.
sol = solveLP(ecModel)

% Growth is reached! Let's make sure we store this functional model
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');

%% STAGE 5: simulation and analysis
% If starting from here, load some basic assets
ecModel = loadEcModel('ecYeastGEMfull.yml');
modelY = loadConventionalGEM();
fluxData = loadFluxData;

% STEP 21 Example of various useful RAVEN functions
% % Set the upper bound of reaction r_0001 to 10.
% ecModel = setParam(ecModel,'ub','r_0001',10);
% % Set the lower bound of reaction r_0001 to 0.
% ecModel = setParam(ecModel,'lb','r_0001',0);
% % Set the objective function to maximize reaction biomassRxn
% ecModel = setParam(ecModel,'obj','r_4041',1);
% % Set the objective function to minimize protein usage
% ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
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

% STEP 22 Selecting objective functions
ecModel = setParam(ecModel,'obj',params.bioRxn,1);
sol = solveLP(ecModel)
disp(['Growth rate reached: ' num2str(abs(sol.f))])
% Set growth lower bound to 99% of the previous value
ecModel = setParam(ecModel,'lb',params.bioRxn,0.99*abs(sol.f));
% Minimize protein pool usage. As protein pool exchange is defined in the
% reverse direction (with negative flux), minimization of protein pool
% usage is computationally represented by maximizing the prot_pool_exchange
% reaction.
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
sol = solveLP(ecModel)
disp(['Minimum protein pool usage: ' num2str(abs(sol.f)) ' ug/gDCW'])

% STEP 23 Compare fluxes from ecModel and starting model
% Constrain with the same conditions to model and ecModel. We now fix the
% observed growth as lower bound ('min' in constrainFluxData) and allow 5%
% variability around the other measured fluxes.
% We know that growth can only reach 0.088, so use this instead of 0.1.
fluxData.grRate(1) = 0.088;
ecModel = constrainFluxData(ecModel,fluxData,1,'min',5);
% Minimize protein pool usage. 
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
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

% STEP 24 Inspect enzyme usage
% Show the result from the earlier simulation, without mapping to
% non-ecModel.
usageData = enzymeUsage(ecModel,solEC.x);
usageReport = reportEnzymeUsage(ecModel,usageData,0.90);

% STEP 25 Perform (ec)FVA
% Perform FVA
[minFluxEc, maxFluxEc] = ecFVA(ecModel, modelY);
[minFluxY, maxFluxY] = ecFVA(modelY, modelY);

% Write results to output file
output = [modelY.rxns, modelY.rxnNames, num2cell([minFluxY, maxFluxY, minFluxEc, maxFluxEc])]';
fID = fopen(fullfile(params.path,'output','ecFVA.tsv'),'w');
fprintf(fID,'%s %s %s %s %s %s\n','rxnIDs', 'rxnNames', 'minFlux', ...
            'maxFlux', 'ec-minFlux', 'ec-maxFlux');
fprintf(fID,'%s %s %g %g %g %g\n',output{:});
fclose(fID);

% Look at flux ranges to indicate reaction-level variability
fluxRange = maxFluxEc - minFluxEc;
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
