% This code requires that Human-GEM is installed

GECKORoot = findGECKOroot();

humanAdapter = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'tutorials', 'tutorial_Human-GEM'));
HumanGEMRoot = HumanGEMAdapter.getHumanGEMRootPath();

ihuman = load(fullfile(HumanGEMRoot, 'model', 'Human-GEM.mat')).ihuman;
% Remove dead-end reactions
modSimp = simplifyModel(ihuman,false,false,true,true);
modExp = modSimp;
[modExp.grRules,skipped] = simplifyGrRules(modExp.grRules,true);

% Full model
%%%%%%%%%%%%

fullECModel = makeEcModel(modExp, false, humanAdapter);
[fullECModel, foundComplex, proposedComplex] = applyComplexData(fullECModel, [], humanAdapter);

% Run DLKcat
fullECModel = findMetSmiles(fullECModel, humanAdapter);
writeDLKcatInput(fullECModel, [], humanAdapter, [], [], true);
runDLKcat(humanAdapter);
kcatListFullDlKcat = readDLKcatOutput(fullECModel, [], humanAdapter);

%Run fuzzy matching
%ecModel         = getECfromDatabase(ecModel, [], 'display', humanAdapter);
ecModel         = getECfromGEM(ecModel);
kcatListFullFuzzy = fuzzyKcatMatching(fullECModel, [], humanAdapter);

mergedKcatListFull = mergeDLKcatAndFuzzyKcats(kcatListFullDlKcat, kcatListFullFuzzy);

fullECModelMerged = selectKcatValue(fullECModel,mergedKcatListFull);
fullECModelMerged = applyKcatConstraints(fullECModelMerged);
%set protein pool constraint
fullECModelMerged = setProtPoolSize(fullECModelMerged, [], [], [], humanAdapter);

%We avoid simulating growth etc. here, it gets complicated with the media

%Light model
%%%%%%%%%%%%

lightECModel = makeEcModel(modExp, true, humanAdapter);
[lightECModel, foundComplex, proposedComplex] = applyComplexData(lightECModel, [], humanAdapter);

%Run DLKcat
lightECModel = findMetSmiles(lightECModel, humanAdapter);
writeDLKcatInput(lightECModel, [], humanAdapter, [], [], true);
runDLKcat(humanAdapter);
kcatListLightDlKcat = readDLKcatOutput(lightECModel, [], humanAdapter);

%Run fuzzy matching
lightECModel = getECfromGEM(lightECModel);
kcatListLightFuzzy = fuzzyKcatMatching(lightECModel, [], humanAdapter);

mergedKcatListLight = mergeDLKcatAndFuzzyKcats(kcatListLightDlKcat, kcatListLightFuzzy);

lightECModelMerged = selectKcatValue(lightECModel, mergedKcatListLight);
lightECModelMerged = applyKcatConstraints(lightECModelMerged);
%set protein pool constraint
lightECModelMerged = setProtPoolSize(lightECModelMerged, [], [], [], humanAdapter);

%We avoid simulating growth etc. here, it gets complicated with the media

%export metabolic genes to figure out the fraction of metabolic enzymes of the total gene expression
%[metabolicGenesHuman,~] = getGenesFromGrRules(ihuman.grRules);
%length(metabolicGenesHuman)%3067
%t = table(metabolicGenesHuman);
%writetable(t, 'metabolicGenesHuman.txt', 'WriteVariableNames', true);















GECKORoot = findGECKOroot();

humanAdapter = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'tutorials', 'tutorial_Human-GEM'));

%test that the path stuff works
%humanAdapterTemp = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'tutorials', 'tutorial_Human-GEMTemp'), false); %should give warning + error
%humanAdapterTemp = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'tutorials', 'tutorial_Human-GEMTemp'), true);
%humanAdapterTemp.getParameters
HumanGEMRoot = HumanGEMAdapter.getHumanGEMRootPath();

ihuman = load(fullfile(HumanGEMRoot, 'model', 'Human-GEM.mat')).ihuman;
modSimp = simplifyModel(ihuman,false,false,true,true); %remove dead-end reactions
modExp = modSimp;
tic
[modExp.grRules,skipped] = simplifyGrRules(modExp.grRules,true);
toc %Elapsed time is 183.114703 seconds.
sum(skipped)%0

tic
fullECModel = makeEcModel(modExp, false, humanAdapter);
toc %Elapsed time is 14.860279 seconds.

tic
lightECModel = makeEcModel(modExp, true, humanAdapter);
toc %Elapsed time is 3.297782 seconds.

lightECModel = getECfromGEM(lightECModel);
%lightECModel2 = getECfromDatabase(lightECModel, [], 'display', humanAdapter);

fullECModel = getECfromGEM(fullECModel);

%Show the differences between the model ec codes and the ones derived from the databases
%sel = ~(cellfun(@isempty, lightECModel.ec.eccodes) & cellfun(@isempty, lightECModel2.ec.eccodes));
%table(lightECModel.ec.eccodes(sel), lightECModel2.ec.eccodes(sel))

tic
%complexInfo = getComplexData('Homo sapiens', humanAdapter);
[lightECModel, foundComplex, proposedComplex] = applyComplexData(lightECModel, [], humanAdapter);
toc %Elapsed time is 62.289315 seconds - a bit slow, should be easy to optimize at some point.
%look at complex II - why does it not match in complexPortal?
%xxx = strcmp(lightECModel.ec.rxns, '001_MAR04652');
%protSel = lightECModel.ec.rxnEnzMat(xxx,:) > 0;
%sum(protSel) %4 - ok
%lightECModel.ec.enzymes(protSel)

tic
%complexInfo = getComplexData('Homo sapiens', humanAdapter);
[fullECModel, foundComplex, proposedComplex] = applyComplexData(fullECModel, [], humanAdapter);
toc %Elapsed time is 62.289315 seconds - a bit slow, should be easy to optimize at some point.


kcatListLightFuzzy = fuzzyKcatMatching(lightECModel, [], humanAdapter);
%lightECModel = selectKcatValue(lightECModel,kcatList);
%lightECModel = applyKcatConstraints(lightECModel);
%constructEquations(lightECModel)
%hmm, look at the prot_pool usage
%find(lightECModel.S(length(lightECModel.mets),:) < -10^12)
%constructEquations(lightECModel, lightECModel.rxns(851))

% Run dlKcat
fullECModel = findMetSmiles(fullECModel, humanAdapter);
testFull = writeDLKcatInput(fullECModel, [], humanAdapter);
runDLKcat(humanAdapter);
kcatListFullDlKcat = readDLKcatOutput(fullECModel, [], humanAdapter);

%fullECModel2 = fullECModel;
%fullECModel2.metSmiles = [];
%fullECModel2.metSmiles = findMetSmiles(fullECModel2.metNames);
%all(strcmp(fullECModel2.metSmiles, fullECModel.metSmiles))
%sel = ~strcmp(fullECModel2.metSmiles, fullECModel.metSmiles)
%fullECModel2.metNames(sel)
%fullECModel.metNames(sel)
%table(fullECModel2.metSmiles(sel),fullECModel.metSmiles(sel))


lightECModel = findMetSmiles(lightECModel, humanAdapter);
testLight = writeDLKcatInput(lightECModel, [], humanAdapter);
runDLKcat(humanAdapter);
kcatListLightDlKcat = readDLKcatOutput(lightECModel, [], humanAdapter);

%now join the fuzzy and dlkcat ckats for light
mergedKcatListLight = mergeDlkcatAndFuzzyKcats(kcatListLightDlKcat, kcatListLightFuzzy);
lightECModelMerged  = selectKcatValue(lightECModel,mergedKcatListLight);
lightECModelMerged = applyKcatConstraints(lightECModelMerged);

% plot the coefficients
coeffs = -full(lightECModelMerged.S(length(lightECModelMerged.mets), 1:(length(lightECModelMerged.rxns)-1))).';
coeffs(coeffs == 0) = [];
length(coeffs)
histogram(log10(coeffs)) % looks like the current coeffs should be multiplied by 1000 - center is around 10^-3

%lightECModel.ec.mw(strcmp(lightECModel.ec.rxns, '001_MAR03875')) %25195, reasonable
%lightECModel.ec.kcat(strcmp(lightECModel.ec.rxns, '001_MAR03875')) %8.8300e-05, very small

% set protein pool constraint
lightECModelMerged = setProtPoolSize(lightECModelMerged, [], humanAdapter);
lightECModelTuned = sensitivityTuning(lightECModel, 0.07, humanAdapter);

% full
kcatListFullFuzzy = fuzzyKcatMatching(fullECModel, [], humanAdapter);
mergedKcatListFull = mergeDlkcatAndFuzzyKcats(kcatListFullDlKcat, kcatListFullFuzzy);

fullECModelMerged = selectKcatValue(fullECModel,mergedKcatListFull);
fullECModelMerged = applyKcatConstraints(fullECModelMerged);

% set protein pool constraint
fullECModelMerged = setProtPoolSize(fullECModelMerged, [], humanAdapter);




% Now merge the kcats from dlkcat and fuzzy
mergedKcatList = mergeDlkcatAndFuzzyKcats(kcatListFullDlKcat, kcatListFullFuzzy);
fullECModelMerged  = selectKcatValue(fullECModelMerged,mergedKcatList);
fullECModelMerged = applyKcatConstraints(fullECModelMerged);

fullECModelMerged = setProtPoolSize(fullECModelMerged, [], humanAdapter);

fullECModelTuned = sensitivityTuning(fullECModelMerged, 0.07, humanAdapter);


constructEquations(fullECModel, 'prot_pool_exchange')
constructEquations(fullECModel, 'MAR12341_EXP_1')
constructEquations(fullECModel, 'draw_prot_Q9Y6K0')

%where is prot_Q9Y6K0 used?
fullECModel.S(strcmp(fullECModel.metNames,'prot_Q9Y6K0'),:)
fullECModel.rxns(3346)

% good example
constructEquations(fullECModel, 'MAR00468_EXP_2')
constructEquations(fullECModel, 'MAR00468_EXP_1')
constructEquations(lightECModel, 'MAR00468')

% look at the coefficients in light and full
sel = startsWith(fullECModel.rxns, 'draw');
constructEquations(fullECModel, fullECModel.rxns(~sel))

kcatList.kcats(strcmp(kcatList.rxns,'001_MAR00468')) %2.7900e-04
kcatListFull.kcats(strcmp(kcatListFull.rxns,'MAR00468_EXP_1'))
kcatListFull.rxns(kcatListFull.kcats ~= 0)

fullECModel.ec.kcat(strcmp(fullECModel.ec.rxns, 'MAR00468_EXP_1'))
lightECModel.ec.kcat(strcmp(lightECModel.ec.rxns, '001_MAR00468')) % ok

