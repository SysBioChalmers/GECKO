%This code requires that Human-GEM is installed

GECKORoot = findGECKOroot();

yeastAdapter = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'userData', 'ecYeastGEM'));

%ihuman = load(fullfile(GECKORoot, 'model', 'Human-GEM.mat')).ihuman;
yeastGEM = load('C:/Code/Components/yeast-GEM/yeast-GEM/model/yeast-GEM.mat').model;

tic
fullECModel = makeEcModel(yeastGEM, false, yeastAdapter);
toc

tic
lightECModel = makeEcModel(yeastGEM, true, yeastAdapter);
toc

lightECModel2 = getECfromDatabase(lightECModel, 'display', [], yeastAdapter); 
lightECModel = getECfromGEM(lightECModel);
fullECModel = getECfromGEM(fullECModel);

%Show the differences between the model ec codes and the ones derived from the databases
%sel = ~(cellfun(@isempty, lightECModel.ec.eccodes) & cellfun(@isempty, lightECModel2.ec.eccodes));
%table(lightECModel.ec.eccodes(sel), lightECModel2.ec.eccodes(sel))

tic
[lightECModel, foundComplex, proposedComplex] = applyComplexData(lightECModel, [], yeastAdapter);
toc 

tic
[fullECModel, foundComplex, proposedComplex] = applyComplexData(fullECModel, [], yeastAdapter);
toc 


kcatListLightFuzzy = fuzzyKcatMatching(lightECModel, [], yeastAdapter);
%lightECModel = selectKcatValue(lightECModel,kcatList);
%lightECModel = applyKcatConstraints(lightECModel);
%constructEquations(lightECModel)
%hmm, look at the prot_pool usage
%find(lightECModel.S(length(lightECModel.mets),:) < -10^12)
%constructEquations(lightECModel, lightECModel.rxns(851))

%run dlKcat
%ModelAdapterManager.setDefaultAdapter(yeastAdapter); %dlkcat functions currently do not well support sending in an adapter - fix later
fullECModel = findMetSmiles(fullECModel, yeastAdapter);
testFull = writeDLKcatInput(fullECModel, [], yeastAdapter);
runDLKcat([], [], yeastAdapter, [], 'C:/Python38/Python38/python.exe', 'C:/Python38/Python38/Scripts/pip.exe');
kcatListFullDlKcat = readDLKcatOutput(fullECModel, [], yeastAdapter);

lightECModel = findMetSmiles(lightECModel, yeastAdapter);
testLight = writeDLKcatInput(lightECModel, [], yeastAdapter);
runDLKcat([], [], yeastAdapter, [], 'C:/Python38/Python38/python.exe', 'C:/Python38/Python38/Scripts/pip.exe');
kcatListLightDlKcat = readDLKcatOutput(lightECModel, [], yeastAdapter);

%now join the fuzzy and dlkcat ckats for light
mergedKcatListLight = mergeDlkcatAndFuzzyKcats(kcatListLightDlKcat, kcatListLightFuzzy);
lightECModelMerged  = selectKcatValue(lightECModel,mergedKcatListLight);
lightECModelMerged = applyKcatConstraints(lightECModelMerged);

%plot the coefficients
coeffs = -full(lightECModelMerged.S(length(lightECModelMerged.mets), 1:(length(lightECModelMerged.rxns)-1))).';
coeffs(coeffs == 0) = [];
length(coeffs)
histogram(log10(coeffs)) %looks like the current coeffs should be multiplied by 1000 - center is around 10^-3


%lightECModel.ec.mw(strcmp(lightECModel.ec.rxns, '001_MAR03875')) %25195, reasonable
%lightECModel.ec.kcat(strcmp(lightECModel.ec.rxns, '001_MAR03875')) %8.8300e-05, very small

%set protein pool constraint
lightECModelMerged = setProtPoolSize(lightECModelMerged, [], yeastAdapter);



lightECModelTuned = sensitivityTuning(lightECModel, 0.07, yeastAdapter);


%full
kcatListFullFuzzy = fuzzyKcatMatching(fullECModel, [], yeastAdapter);

mergedKcatListFull = mergeDlkcatAndFuzzyKcats(kcatListFullDlKcat, kcatListFullFuzzy);

fullECModelMerged = selectKcatValue(fullECModel,mergedKcatListFull);
fullECModelMerged = applyKcatConstraints(fullECModelMerged);
%set protein pool constraint
fullECModelMerged = setProtPoolSize(fullECModelMerged, [], yeastAdapter);
tunedFullModel = BayesianSensitivityTuning(fullECModelMerged,yeastAdapter,150);

tunedLightModel = BayesianSensitivityTuning(lightECModelMerged,yeastAdapter,150);

%fullECModelTuned = sensitivityTuning(fullECModelMerged, 0.07, yeastAdapter);




constructEquations(fullECModel, 'prot_pool_exchange')
constructEquations(fullECModel, 'MAR12341_EXP_1')
constructEquations(fullECModel, 'draw_prot_Q9Y6K0')
%where is prot_Q9Y6K0 used?
fullECModel.S(strcmp(fullECModel.metNames,'prot_Q9Y6K0'),:)
fullECModel.rxns(3346)

%good example
constructEquations(fullECModel, 'MAR00468_EXP_2')
constructEquations(fullECModel, 'MAR00468_EXP_1')
constructEquations(lightECModel, 'MAR00468')

%look at the coefficients in light and full
sel = startsWith(fullECModel.rxns, 'draw');
constructEquations(fullECModel, fullECModel.rxns(~sel))

kcatList.kcats(strcmp(kcatList.rxns,'001_MAR00468')) %2.7900e-04
kcatListFull.kcats(strcmp(kcatListFull.rxns,'MAR00468_EXP_1'))
kcatListFull.rxns(kcatListFull.kcats ~= 0)

fullECModel.ec.kcat(strcmp(fullECModel.ec.rxns, 'MAR00468_EXP_1'))
lightECModel.ec.kcat(strcmp(lightECModel.ec.rxns, '001_MAR00468'))%ok

