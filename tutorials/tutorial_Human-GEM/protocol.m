GECKOInstaller.install
checkInstallation

%% STAGE 1: expansion from a starting metabolic model to an ecModel structure

% STEP 1 Set modelAdapter
geckoRoot = findGECKOroot;
ModelAdapterManager.setDefaultAdapterFromPath(fullfile(geckoRoot,'tutorials','ecYeastGEM')); 
ModelAdapter = ModelAdapterManager.getDefaultAdapter();
params = ModelAdapter.getParameters();

% STEP 2 Load conventional GEM
modelY = loadConventionalGEM();

% STEP 3 Prepare ecModel
[ecModelFull, noUniprot] = makeEcModel(modelY,false);
[ecModelLight, noUniprot] = makeEcModel(modelY,true);
ecModel = ecModelFull;
doc makeEcModel

% STEP 4 Annotate with complex data
complexInfo = getComplexData(); % No need to run, as data is already in adapter folder
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

% STEP 5 Store model in YAML format
saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');
ecModel=loadEcModel('ecYeastGEMfull.yml');

%% STAGE 2: integration of kcat into the ecModel structure

% STEP 6 Fuzzy matching with BRENDA
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel,noEC);
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% STEP 7 DLKcat prediction through machine learning
runDLKcat();
kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 8 Combine kcat from BRENDA and DLKcat
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

% STEP 9 Take kcatList and populate edModel.ec.kcat
ecModel  = selectKcatValue(ecModel, kcatList_merged);

% STEP 10 Apply custom kcat values
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);
ecModel  = applyKcatConstraints(ecModel);

% STEP 11 Get kcat values across isoenzymes
ecModel = getKcatAcrossIsoenzymes(ecModel);

% STEP 12 Get standard kcat
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
ecModel = applyKcatConstraints(ecModel);

% STEP 14 Set upper bound of protein pool
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;

ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

%% STAGE 3: model tuning
ecModel = setParam(ecModel,'lb','r_1714',-1000);
sol = solveLP(ecModel,1)
printFluxes(ecModel, sol.x)

% STEP 15 Relax protein pool constraint
protPoolIdx = strcmp(ecModel.rxns, 'prot_pool_exchange');
ecModel.lb(protPoolIdx) = -1000;
sol = solveLP(ecModel,1)
ecModel.lb(protPoolIdx) = sol.x(protPoolIdx);

% STEP 16 Sensitivity tuning
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

%% STAGE 4 integration of proteomics data into the ecModel
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

% STEP 20 Protein concentrations are flexibilized (increased), until the
[ecModel, flexProt] = flexibilizeProtConcs(ecModel,0.1,10);

modelY = constrainFluxData(modelY,fluxData);
sol = solveLP(modelY)

sol = solveLP(ecModel)

saveEcModel(ecModel,ModelAdapter,'yml','ecYeastGEMfull');

%% STAGE 5: simulation and analysis
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
ecModel = setParam(ecModel,'lb',params.bioRxn,0.99*abs(sol.f));
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
sol = solveLP(ecModel)
disp(['Minimum protein pool usage: ' num2str(abs(sol.f)) ' mg/gDCW'])

% STEP 23 Compare fluxes from ecModel and starting model
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

% STEP 24 Inspect enzyme usage
usageData = enzymeUsage(ecModel,solEC.x);
usageReport = reportEnzymeUsage(ecModel,usageData,0.90);

% STEP 25 Perform (ec)FVA
[minFluxEc, maxFluxEc] = ecFVA(ecModel, modelY);
[minFluxY, maxFluxY] = ecFVA(modelY, modelY);

% Write results to output file
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
