
GECKOInstaller.install
checkInstallation

%% STAGE 1: expansion from a starting metabolic model to an ecModel structure

disp('STEP 1 Set modelAdapter');
geckoRoot = findGECKOroot;
ModelAdapterManager.setDefaultAdapterFromPath(fullfile(geckoRoot,'tutorials', 'tutorial_yeast-GEM', 'YeastGEMAdapter.m')); 
ModelAdapter = ModelAdapterManager.getDefaultAdapter();
params = ModelAdapter.getParameters();

disp('STEP 2 Load conventional GEM');
modelY = loadConventionalGEM();

disp('STEP 3 Prepare ecModel');
[ecModelFull, noUniprot] = makeEcModel(modelY,false);
[ecModelLight, noUniprot] = makeEcModel(modelY,true);
ecModel = ecModelFull;
doc makeEcModel

disp('STEP 4 Annotate with complex data');
complexInfo = getComplexData(); % No need to run, as data is already in adapter folder
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

disp('STEP 5 Store model in YAML format');
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');
ecModel=loadEcModel('ecYeastGEMfull.yml');

%% STAGE 2: integration of kcat into the ecModel structure

disp('STEP 6 Fuzzy matching with BRENDA');
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel,noEC);
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

disp('STEP 7 DLKcat prediction through machine learning');
runDLKcat();
kcatList_DLKcat = readDLKcatOutput(ecModel);

disp('STEP 8 Combine kcat from BRENDA and DLKcat');
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

disp('STEP 9 Take kcatList and populate edModel.ec.kcat');
ecModel  = selectKcatValue(ecModel, kcatList_merged);

disp('STEP 10 Apply custom kcat values');
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);
ecModel  = applyKcatConstraints(ecModel);

disp('STEP 11 Get kcat values across isoenzymes');
ecModel = getKcatAcrossIsoenzymes(ecModel);

disp('STEP 12 Get standard kcat');
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

disp('STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S');
ecModel = applyKcatConstraints(ecModel);

disp('STEP 14 Set upper bound of protein pool');
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;

ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

%% STAGE 3: model tuning
ecModel = setParam(ecModel,'lb','r_1714',-1000);
sol = solveLP(ecModel,1)
printFluxes(ecModel, sol.x)

disp('STEP 15 Relax protein pool constraint');
protPoolIdx = strcmp(ecModel.rxns, 'prot_pool_exchange');
ecModel.lb(protPoolIdx) = -1000;
sol = solveLP(ecModel,1)
ecModel.lb(protPoolIdx) = sol.x(protPoolIdx);

disp('STEP 16 Sensitivity tuning');
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

%% STAGE 4 integration of proteomics data into the ecModel
protData = loadProtData(3); %Number of replicates
ecModel = fillProtConcs(ecModel,protData);
ecModel = constrainProtConcs(ecModel);

% STEP 18 Update protein pool
% The protein pool reaction will be constraint by the remaining, unmeasured
% enzyme content. This is calculated by subtracting the sum of 
% ecModel.ec.concs from the condition-specific total protein content. The
% latter is stored together with the flux data that otherwise will be used
% in Step 19.
fluxData = loadFluxData();
ecModel = updateProtPool(ecModel,fluxData.Ptot(1));

% STEP 19 Load flux data
% Matching the proteomics sample(s), condition-specific flux data needs to
% be loaded to constrain the model. This was already loaded in Step 18 for
% gathering Ptot, but is repeated here nonetheless. Flux data is read from
% /data/fluxData.tsv.
fluxData = loadFluxData();
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose'); % Use first condition
sol = solveLP(ecModel); % To observe if growth was reached
disp(['Growth rate that is reached: ' num2str(abs(sol.f))])

% STEP 20 Protein concentrations are flexibilized (increased), until the
% intended growth rate is reached. This is condition-specific, so the
% intended growth rate is gathered from the fluxData structure.
[ecModel, flexProt] = flexibilizeProtConcs(ecModel,fluxData.grRate(1),10);

modelY = constrainFluxData(modelY,fluxData);
sol = solveLP(modelY)

sol = solveLP(ecModel)

saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');

%% STAGE 5: simulation and analysis
ecModel = loadEcModel('ecYeastGEMfull.yml');
modelY = loadConventionalGEM();
fluxData = loadFluxData;

disp('STEP 21 Example of various useful RAVEN functions');
% skipped here

disp('STEP 22 Selecting objective functions');
ecModel = setParam(ecModel,'obj',params.bioRxn,1);
sol = solveLP(ecModel)
disp(['Growth rate reached: ' num2str(abs(sol.f))])
ecModel = setParam(ecModel,'lb',params.bioRxn,0.99*abs(sol.f));
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
sol = solveLP(ecModel)
disp(['Minimum protein pool usage: ' num2str(abs(sol.f)) ' mg/gDCW'])

disp('STEP 23 Compare fluxes from ecModel and starting model');
fluxData.grRate(1) = 0.088;
ecModel = constrainFluxData(ecModel,fluxData,1,'min',5);
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
solEC = solveLP(ecModel,1)

modelY = constrainFluxData(modelY,fluxData,1,'min',5);
sol = solveLP(modelY,1)

[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, modelY, solEC.x);

printFluxes(modelY,[sol.x mappedFlux])
printFluxes(modelY,[sol.x mappedFlux],false)
ratioFlux = mappedFlux ./ sol.x;
ratioFlux(isnan(ratioFlux)) = 0; % Divisions by zero give NaN, reset to zero.
printFluxes(modelY,ratioFlux,false)

disp('STEP 24 Inspect enzyme usage');
usageData = enzymeUsage(ecModel,solEC.x);
usageReport = reportEnzymeUsage(ecModel,usageData,0.90);

disp('STEP 25 Perform (ec)FVA');
[minFluxEc, maxFluxEc] = ecFVA(ecModel, modelY);
[minFluxY, maxFluxY] = ecFVA(modelY, modelY);

output = [modelY.rxns, modelY.rxnNames, num2cell([minFluxY, maxFluxY, minFluxEc, maxFluxEc])]';
fID = fopen(fullfile(params.path,'output','ecFVA.tsv'),'w');
fprintf(fID,'%s %s %s %s %s %s\n','rxnIDs', 'rxnNames', 'minFlux', ...
            'maxFlux', 'ec-minFlux', 'ec-maxFlux');
fprintf(fID,'%s %s %g %g %g %g\n',output{:});
fclose(fID);

fluxRange = maxFluxEc - minFluxEc;
fluxRangeY = maxFluxY - minFluxY;

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
