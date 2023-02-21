%this is tcm0001

GECKORoot = findGECKOroot();

yeastAdapter = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'userData', 'ecYeastGEM'));
yeastGEM = importModel(fullfile(yeastAdapter.getParameters().path,'models','yeast-GEM.xml'));

%Full model
%%%%%%%%%%%%

fullECModel = makeEcModel(yeastGEM, false, yeastAdapter);
[fullECModel, foundComplex, proposedComplex] = applyComplexData(fullECModel, [], yeastAdapter);

%Run DLKcat
fullECModel = findMetSmiles(fullECModel, yeastAdapter);
writeDLKcatInput(fullECModel, [], yeastAdapter, [], [], true);
runDLKcat([], yeastAdapter);
kcatListFullDlKcat = readDLKcatOutput(fullECModel, [], yeastAdapter);

%Run fuzzy matching
%fullECModel = getECfromDatabase(fullECModel, 'display', [], yeastAdapter);
fullECModel = getECfromGEM(fullECModel);
kcatListFullFuzzy = fuzzyKcatMatching(fullECModel, [], yeastAdapter);

mergedKcatListFull = mergeDLKcatAndFuzzyKcats(kcatListFullDlKcat, kcatListFullFuzzy);

fullECModelMerged = selectKcatValue(fullECModel,mergedKcatListFull);
fullECModelMerged = applyKcatConstraints(fullECModelMerged);
%set protein pool constraint
fullECModelMerged = setProtPoolSize(fullECModelMerged, [], [], [], yeastAdapter);

%growth before tuning
%we need to increase the lb and ub to reach max growth
fullECModelMerged.lb(fullECModelMerged.lb == -1) = -1000;
fullECModelMerged.lb(fullECModelMerged.lb == -10) = -1000;
fullECModelMerged.ub(fullECModelMerged.ub == 1) = 1000;
sol = solveLP(fullECModelMerged, 1);
growthBef = -sol.f;
fullECModelTuned = sensitivityTuning(fullECModelMerged, 0.4, yeastAdapter);
sol = solveLP(fullECModelTuned, 1);
growthAfter = -sol.f;
disp(['Growth before: ' num2str(growthBef) ', Growth after: ' num2str(growthAfter)])

%Light model
%%%%%%%%%%%%

lightECModel = makeEcModel(yeastGEM, true, yeastAdapter);
[lightECModel, foundComplex, proposedComplex] = applyComplexData(lightECModel, [], yeastAdapter);

%Run DLKcat
lightECModel = findMetSmiles(lightECModel, yeastAdapter);
writeDLKcatInput(lightECModel, [], yeastAdapter, [], [], true);
runDLKcat([], yeastAdapter);
kcatListLightDlKcat = readDLKcatOutput(lightECModel, [], yeastAdapter);

%Run fuzzy matching
lightECModel = getECfromGEM(lightECModel);
kcatListLightFuzzy = fuzzyKcatMatching(lightECModel, [], yeastAdapter);

mergedKcatListLight = mergeDLKcatAndFuzzyKcats(kcatListLightDlKcat, kcatListLightFuzzy);

lightECModelMerged = selectKcatValue(lightECModel, mergedKcatListLight);
lightECModelMerged = applyKcatConstraints(lightECModelMerged);
%set protein pool constraint
lightECModelMerged = setProtPoolSize(lightECModelMerged, [], [], [], yeastAdapter);

%growth before tuning
%we need to increase the lb and ub to reach max growth
lightECModelMerged.lb(lightECModelMerged.lb == -1) = -1000;
lightECModelMerged.lb(lightECModelMerged.lb == -10) = -1000;
lightECModelMerged.ub(lightECModelMerged.ub == 1) = 1000;
sol = solveLP(lightECModelMerged, 1);
growthBef = -sol.f;
lightECModelTuned = sensitivityTuning(lightECModelMerged, 0.4, yeastAdapter);
sol = solveLP(lightECModelTuned, 1);
growthAfter = -sol.f;
disp(['Growth before: ' num2str(growthBef) ', Growth after: ' num2str(growthAfter)])

