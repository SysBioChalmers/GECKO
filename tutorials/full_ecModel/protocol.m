% This file accompanies the GECKO3 Nature Protocols paper (DOI TO BE ADDED).
%
% The function of this script is to demonstrate the reconstruction and
% analysis of a *full* ecModel. As example, it here uses the yeast-GEM
% model of Saccharomyces cerevisiae as starting point. However, this script
% does not claim to construct a "production-ready" ecYeastGEM model:
% dependent on how you intend to use the ecModel it may require additional
% curation and evaluation.
%
% DO NOT USE THE ECMODEL GENERATED HERE OUTSIDE OF THIS TUTORIAL.
%
% In comparison to the published GECKO3 Nature Protocols paper, this script
% might have more up-to-date descriptions about the capabilities and
% functions, which were introduced after the Nature Protocols paper was
% published. This script will have additional commands and analyses that
% not included as such in the paper. This script should always work with
% the most recent GECKO3 release.

%% Preparation stage for ecModel reconstruction
% - Install GECKO & RAVEN
%   - GECKO can be installed via cloning or direct download of ZIP file.
%     See installation instructions in the README.md:
%     https://github.com/SysBioChalmers/GECKO/tree/main#installation
%   - Add the appropriate GECKO (sub)folders to MATLAB path:
GECKOInstaller.install

%   - Install RAVEN by following the installation instructions:
%     https://github.com/SysBioChalmers/RAVEN/wiki/Installation
%   - The installation of Gurobi as LP solver is highly recommended
checkInstallation % Confirm that RAVEN is functional, should be 2.7.12 or later.

% - Initiate a basic structure of files and folders for your intended
%   project. This includes a copy of the template adapter.
%   The next line is commented out as the project structure is already
%   available in GECKO/tutorials/full_ecModel.
%startGECKOproject()

% - Find a high-quality GEM of your species of interest. ecYeastGEM is
%   based on yeast-GEM https://github.com/SysBioChalmers/yeast-GEM/releases
% - Release v8.6.2 of yeast-GEM is also distributed with GECKO at 
%   userData/ecYeastGEM/model/yeast-GEM.xml
% - Modify the model adapter (at userData/ecYeastGEM/ecYeastGEMadapter.m)
%   to contain organism- and model-specific parameters.

%% STAGE 1: expansion from a starting metabolic model to an ecModel structure
% STEP 1 Set modelAdapter
adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel','YeastGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);

% With the above line, we set the YeastGEMAdapter as the default adapter
% from here onwards, which means that any GECKO function that requires a
% model adapter as input will use YeastGEMAdapter. This is until you have
% set a new default adapter, or closed MATLAB.

% We can also load the ModelAdapter into the MATLAB Workspace, so that you
% can quickly see its content, for instance the parameters.
ModelAdapter = ModelAdapterManager.getDefault();
params = ModelAdapter.getParameters();
% However, changes should not be made in this structure. Making a change in
% ModelAdapter in the Workspace has no effect (unless the ModelAdapter
% is explicitly provided as a function input), rather you should make the
% change in the adapter file and run anew ModelAdapterManager.setDefault().

% STEP 2 Load conventional yeast-GEM
% If the location to the conventional GEM was already set in the modelAdapter,
% as obj.param.convGEM, then loadConventionalGEM can directly be used. If
% the model is stored somewhere else (and not specified in obj.param.convGEM),
% you can also use RAVEN's importModel(). In that case you will never use
% loadConventionalGEM and the obj.param.convGEM never has to be specified.
model = loadConventionalGEM();
% model = importModel(fullfile(geckoRoot,'tutorials','full_ecModel','models','yeast-GEM.xml'));

% STEP 3 Prepare ecModel
% We will make a full GECKO ecModel. For an example of reconstructing a 
% light GECKO ecModel, see tutorials/light_ecModel
[ecModel, noUniprot] = makeEcModel(model,false);
% Note that noUniprot is empty: for all genes a match could be find in the
% Uniprot dataset

% It can be helpful to read the function documentation: to read its
% function and what is the order of inputs and outputs.
doc makeEcModel

% STEP 4 Annotate with complex data
% As Saccharomyces cerevisiae is available on Complex Portal
% (https://www.ebi.ac.uk/complexportal/), its data is included in
% ecYeastGEM.
% The next line is commented out, as the data is already available in the
% full_ecModel/data folder.
%complexInfo = getComplexData();
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

% complexInfo can be given as second input, but not needed here, as it will
% read the file that was written by getComplexData

% STEP 5 Store model in YAML format
% In this script there is code at the end of each stage to store a copy of
% the ecModel. However, only the stage 3 ecModel is kept and distributed
% together with GECKO (as is shown at the end of stage 3).
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEM_stage1');

%% STAGE 2: integration of kcat into the ecModel structure
%ecModel=loadEcModel('ecYeastGEM_stage1.yml'); % Uncomment if you want to reload model

% Decide which kcat source to use. In the steps below, all options are
% shown. Not all are required, it is up to the user to decide which ones
% they want to use.

% STEP 6 Fuzzy matching with BRENDA
% Requires EC numbers, which are here first taken from the starting model,
% with the missing ones taken from Uniprot & KEGG databases.
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel,noEC);
% However, the EC numbers in the yeast-GEM model have not been thoroughly
% curated. Instead, we will take database-derived EC numbers for all
% reactions.
ecModel         = getECfromDatabase(ecModel);

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
% applyKcatConstraints should be run. This is done in STEP 13.

% STEP 11 Get kcat values across isoenzymes
ecModel = getKcatAcrossIsoenzymes(ecModel);

% STEP 12 Get standard kcat
% Assign a protein cost to reactions without gene assocation. These
% reactions are identified as those with empty entry in ecModel.grRules.
% The following reactions are exempted:
% A Exchange reactions: exchanging a metabolite across the model boundary,
%   not representing a real enzymatic reaction.
% B Transport reactions: transporting a metabolite with the same name from
%   one compartment to another. Real transport reactions should already be
%   annotated with grRules, so that the remaining non-annotated reactions
%   are mostly representing diffusion or pseudotransport processes such as
%   vesicles moving from ER to Golgi. While proteins are involved in such
%   processes, they are not catalyzed by enzymes.
% C Pseudoreactions: any other reaction that should not be considered to be
%   catalyzed by an enzyme. getStandardKcat recognizes these from the
%   reaction name contaning "pseudoreaction".
% D Custom list of non-enzyme reactions: if the above approaches does not
%   correctly identify all non-enzyme reactions that should be ignored by
%   getStandardKcat, /data/pseudoRxns.tsv can be specified in the adapter
%   folder, containing the relevant reaction identifiers.
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
% While STEP 9-12 add or modify values in ecModel.ec.kcat, these contraints
% are not yet applied to the S-matrix: the enzyme is not yet set as pseudo-
% substrate with the -MW/kcat stoichiometry. Now, applyKcatConstraints will
% take the values from ecModel.ec.kcat, considers any complex/subunit data
% that is tracked in ecModel.ec.rxnEnzMat, together with the MW in
% ecModel.ec.mw, and uses this to modify the enzyme usage coefficients
% directly in ecModel.S. Any time a change is made to the .kcat, .rxnEnzMat
% or .mw fields, the applyKcatConstraints function should be run again to
% reapply the new constraints onto the metabolic model.
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
% is available via PAXdb (https://pax-db.org/).
% calculateFfactor can be used to estimate the f-factor.
%f = calculateFfactor(ecModel); % Optional
ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

% Note that at a later stage (after stage 3), the sigma factor be further
% adjusted with sigmaFitter, to get a model that is able to reach a
% particular maximum growth rate. This will not be done here, as we first
% need to tune the kcat values in Stage 3.

saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEM_stage2');

%% STAGE 3: model tuning
%ecModel=loadEcModel('ecYeastGEM_stage2.yml'); % Uncomment if you want to reload model

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

% Inspect the tunedKcats structure in table format
struct2table(tunedKcats)

% STEP 17 Curate kcat values based on kcat tuning
% As example, the kcat of 5'-phosphoribosylformyl glycinamidine synthetase
% (reaction r_0079) was increased from 0.05 to 5. Inspecting the kcat
% source might help to determine if this is reasonable. 
rxnIdx = find(strcmp(kcatList_merged.rxns,'r_0079'));
doc fuzzyKcatMatching % To check the meaning of wildcardLvl and origin
kcatList_merged.wildcardLvl(rxnIdx) % 0: no EC number wildcard
kcatList_merged.origin(rxnIdx) % 4: any organism, any substrate, kcat
kcatList_merged.eccodes(rxnIdx) % EC number 6.3.5.3

% On BRENDA https://www.brenda-enzymes.org/enzyme.php?ecno=6.3.5.3#TURNOVER%20NUMBER%20[1/s]
% The kcat value is from E. coli with NH4+ as substrate. The reaction
% normally uses glutamine, so this kcat value might be misleading.
% Inspecting the abstract of the paper that is reporting this value https://pubmed.ncbi.nlm.nih.gov/2659070/
% actually states "and NH3 can replace glutamine as a nitrogen donor with a
% Km = 1 M and a turnover of 3 min-1 (2% glutamine turnover)". The paper
% also reports a specific activity that can be used instead:
% https://www.brenda-enzymes.org/enzyme.php?ecno=6.3.5.3#SPECIFIC%20ACTIVITY%20[%C2%B5mol/min/mg]

% Convert specific activity of 2.15 umol/min/mg protein, where the protein
% has a molecular weight of 148905 Da, to kcat in /sec:

enzMW = ecModel.ec.mw(strcmp(ecModel.ec.enzymes,'P38972')); % Get MW of the enzyme
convKcat = 2.15; % umol/min/mg protein, same as mmol/min/g protein
convKcat = convKcat / 1000; % mol/min/g protein
convKcat = convKcat / 60; % mol/sec/g protein
convKcat = convKcat * enzMW % mol/sec/mol protein, same as 1/sec.

% New kcat is 5.3358, which is not far away from the tuned kcat of 5.

% Another kcat that drastically changed is tryptophan synthase (r_1055),
% with kcat increased from 0.023 to 2.3. Again, check what the source was:
rxnIdx = find(strcmp(kcatList_merged.rxns,'r_1055'));
kcatList_merged.wildcardLvl(rxnIdx) % 0: no EC number wildcard
kcatList_merged.eccodes(rxnIdx) % EC number 4.2.1.20
kcatList_merged.origin(rxnIdx) % 2: any organism, correct substrate, kcat
% Searching for the EC number in GECKO/databases/max_KCAT.txt, it appears
% that the kcat value from Zea mays was taken, as this was the most closely
% related organism with reported kcat. However, corn is not that close to
% yeast, while at the same time there is a S. cerevisiae specific activity
% reported on BRENDA: https://www.brenda-enzymes.org/enzyme.php?ecno=4.2.1.20&organism%5B%5D=Saccharomyces+cerevisiae#SPECIFIC%20ACTIVITY%20[%C2%B5mol/min/mg]

% Again convert this specific activity of 1.69 umol/min/mg protein, with
% protein weight of 76626 Da, to kcat in /sec:

convKcat = 1.69; % umol/min/mg protein, same as mmol/min/g protein
convKcat = convKcat / 1000; % mol/min/g protein
convKcat = convKcat / 60; % mol/sec/g protein
convKcat = convKcat * 76626 % mol/sec/mol protein, same as 1/sec.

% New kcat is 2.16, which is not too far from what sensitivity tuning
% suggested at 2.23.

saveEcModel(ecModel,[],'yml','ecYeastGEM_stage3');

% This functional ecModel will also be kept in the GECKO GitHub
saveEcModel(ecModel,[],'yml','ecYeastGEM');

%% STAGE 4 integration of proteomics data into the ecModel
%ecModel=loadEcModel('ecYeastGEM_stage3.yml'); % Uncomment if you want to reload model

% STEP 18 Load proteomics data and constrain ecModel 
protData = loadProtData(3); %Number of replicates
ecModel = fillProtConcs(ecModel,protData);
ecModel = constrainProtConcs(ecModel);

% STEP 19 Update protein pool
% The protein pool reaction will be constraint by the remaining, unmeasured
% enzyme content. This is calculated by subtracting the sum of 
% ecModel.ec.concs from the condition-specific total protein content. The
% latter is stored together with the flux data that otherwise will be used
% in Step 20.
fluxData = loadFluxData();
ecModel = updateProtPool(ecModel,fluxData.Ptot(1));

% STEP 20 Load flux data
% Matching the proteomics sample(s), condition-specific flux data needs to
% be loaded to constrain the model. This was already loaded in Step 18 for
% gathering Ptot, but is repeated here nonetheless. Flux data is read from
% /data/fluxData.tsv.
fluxData = loadFluxData();
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose'); % Use first condition
sol = solveLP(ecModel); % To observe if growth was reached
disp(['Growth rate that is reached: ' num2str(abs(sol.f))])
% Growth rate of 0.1 is by far not reached, flexibilize protein
% concentrations

% STEP 21 Protein concentrations are flexibilized (increased), until the
% intended growth rate is reached. This is condition-specific, so the
% intended growth rate is gathered from the fluxData structure.
[ecModel, flexProt] = flexibilizeProtConcs(ecModel,fluxData.grRate(1),10);

% Neither individual protein levels nor total protein pool are limiting
% growth. Test whether the starting model is able to reach 0.1.
model = constrainFluxData(model,fluxData);
sol = solveLP(model)

% It also only reaches 0.0889! So the metabolic network would not be able
% to adhere to all measured constraints. Perhaps there is something
% incorrect with the measurements? Regardless, the ecModel is now able to
% reach about 0.0889, which will be fine for now.
sol = solveLP(ecModel)

% Inspect the flexibilized proteins
struct2table(flexProt)
% Inspect which proteins have theri amount changed the most
protChange = flexProt.flexConcs./flexProt.oldConcs;
[ratio, index] = max(protChange)
flexProt.uniprotIDs(index)

% Growth is reached! Let's make sure we store this functional model
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEM_stage4');

%% STAGE 5: simulation and analysis
% STEP 22 Example of various useful RAVEN functions
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

% STEP 23 Simulate Crabtree effect with protein pool

% (Re)load the ecModel without proteomics integration
ecModel = loadEcModel('ecYeastGEM.yml');
% We will soon run a custom plotCrabtree function that is kept in the code
% subfolder. To run this function we will need to navigate into the folder
% where it is stored, but we will navigate back to the current folder
% afterwards.
currentFolder = pwd;
cd(fullfile(ModelAdapter.params.path,'code'))
[fluxes, gRate] = plotCrabtree(ecModel);
% fluxes has all the predicted fluxes, while gRate is a vector with the
% corresponding growth rates that were simulated, as visualized on the
% x-axis in the graph.
% The plot will also be saved in the output subfolder.
saveas(gcf,fullfile(ModelAdapter.params.path,'output','crabtree.tiff'))

% The two graphs show (left:) exchange fluxes from simulations (lines) and
% experiments (circles, from doi:10.1128/AEM.64.11.4226-4233.1998); and
% (right:) the fraction of protein pool that is used by enzymes. In the
% left graph, the y-axis indicates absolute fluxes, so that glucose uptake
% and CO2 excretion both have positive numbers. The model simulation
% demonstrates the Crabtree-effect: at increasing growth rates yeast
% switches from respiration to fermentation, and this occurs when the
% protein pool becomes fully used and thereby limiting. The shift away from
% respiration is most clearly shown by reduced oxygen uptake and increased
% ethanol excretion.

% For comparison, make a similar Crabtree plot for a conventional GEM
% Set protein pool to infinite, to mimic a conventional GEM
ecModel_infProt=setProtPoolSize(ecModel,Inf);
plotCrabtree(ecModel_infProt);
saveas(gcf,fullfile(ModelAdapter.params.path,'output','crabtree_infProt.tiff'))
% It is obvious that no total protein constraint is reached, and Crabtree
% effect is not observed.

% Perform the Crabtree simulation on the pre-Step 16 ecModel (where kcat
% sensitivity tuning has not yet been applied).
ecModel_preTuning = loadEcModel('ecYeastGEM_stage2.yml');
ecModel_preTuning = setParam(ecModel_preTuning,'lb','r_1714',-1000);
plotCrabtree(ecModel_preTuning);
saveas(gcf,fullfile(ModelAdapter.params.path,'output','crabtree_preStep16.tiff'))
% Without kcat tuning, the model gets constrained too early (at too low
% growth rates), which means that no solutions exist at high growth rates.

% Going back to the Crabtree-effect in the post-sensitivity tuning model.
% While the plot shows the same overall trend as the experimental data,
% the protein pool seems to become limiting too early. The slope of the
% predicted glucose uptake rate starts to go up after a growth rate of 0.2,
% while in the experimental data it seems like this should happen around
% 0.28 instead. To see why the protein pool is becoming limiting too early,
% we can look at the most limiting enzyme at growth rate 0.25.

% At which position in the fluxes vector is growth rate 0.25?
find(round(gRate,3) == 0.25)
% Gather enzyme usage data at growth rate 0.25
usageData = enzymeUsage(ecModel, fluxes(:,11));
% Prepare usage report
usageReport = reportEnzymeUsage(ecModel,usageData);
% Inspect the topAbsUsage table, with the top-10 absolute usages. The
% highCapUsage table is empty, as no proteomics is integrated and no
% capacity usages can be calculated.

usageReport.topAbsUsage
% The top used enzyme is the ADP/ATP transporter, using 12% of the total
% available protein pool. The kcat of this reaction was estimated by the
% standardKcat function, and assigned a value of 11. This is most likely a
% too-low value: resulting in too-high enzyme demand. Instead of the
% standard kcat, we can look into literature to find more realistic
% estimates. In https://pubs.acs.org/doi/full/10.1021/bi960668j, Figure 3A,
% a maximum transport activity of 8.3 mmol/min/g protein is reported, while
% elsewhere in the article is noted that the purity of the enzyme
% preparation was around 50%. Conversion from specific activity to kcat
% gives: 9.2166 /sec.
AAt = 8.3*2;        % mmol/min/g protein
AAt = AAt / 1000;   % mol/min/g protein
AAt = AAt / 60;     % mol/sec/g protein
AAt = AAt * 33313;  % 1/sec (as MW of the transporter is 33313 Da [g/mol])

% Set this kcat value for the reaction in both directions
ecModel = setKcatForReactions(ecModel,'r_1110',AAt);
ecModel = setKcatForReactions(ecModel,'r_1110_REV',AAt);
ecModel = applyKcatConstraints(ecModel);

% Inspect the new Crabtree plot: not much difference
plotCrabtree(ecModel);

% Instead, fit the sigma factor
[ecModel,sigma,RMSE] = sigmaFitterCrabtree(ecModel);
fprintf('New sigma factor: %.2f\n', sigma)
plotCrabtree(ecModel);

% Inspect maximum growth rate
ecModel = setParam(ecModel,'obj','r_2111',1);
sol=solveLP(ecModel);
disp(['Growth rate reached: ' num2str(abs(sol.f))])

% If starting from here, load some basic assets
model = loadConventionalGEM();
fluxData = loadFluxData;

%% === FROM HERE NOT FULLY TESTED/UPDATED YET ===

% STEP 24 Selecting objective functions
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
disp(['Minimum protein pool usage: ' num2str(abs(sol.f)) ' mg/gDCW'])

% STEP 25 Compare fluxes from ecModel and starting model
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
model = constrainFluxData(model,fluxData,1,'min',5);
sol = solveLP(model,1)
% Note that with the 5% of variability, the model can reach higher growth
% rate than 0.0885, as it is allowed to take up a little bit more glucose
% this time.

% Map the ecModel fluxes back to the non-ecModel
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, model, solEC.x);

% Print the fluxes next to each other, showing only exchange reactions
printFluxes(model,[sol.x mappedFlux])
% Print the fluxes next to each other, showing all reactions
printFluxes(model,[sol.x mappedFlux],false)
% Look at ratio of ecModel / non-ecModel
ratioFlux = mappedFlux ./ sol.x;
ratioFlux(isnan(ratioFlux)) = 0; % Divisions by zero give NaN, reset to zero.
printFluxes(model,ratioFlux,false)

% STEP 26 Inspect enzyme usage
% Show the result from the earlier simulation, without mapping to
% non-ecModel.
usageData = enzymeUsage(ecModel,solEC.x);
usageReport = reportEnzymeUsage(ecModel,usageData,0.90);

% STEP 25 Perform (ec)FVA
% Perform FVA
[minFluxEc, maxFluxEc] = ecFVA(ecModel, model);
[minFluxY, maxFluxY] = ecFVA(model, model);

% Write results to output file
output = [model.rxns, model.rxnNames, num2cell([minFluxY, maxFluxY, minFluxEc, maxFluxEc])]';
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
