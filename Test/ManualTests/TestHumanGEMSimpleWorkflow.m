%This code requires that Human-GEM is installed

GECKORoot = findGECKOroot();

humanAdapter = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'userData', 'ecHumanGEM'));

%test that the path stuff works
%humanAdapterTemp = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'userData', 'ecHumanGEMTemp'), false); %should give warning + error
%humanAdapterTemp = ModelAdapterManager.getAdapterFromPath(fullfile(GECKORoot, 'userData', 'ecHumanGEMTemp'), true);
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
lightECModel2 = getECfromDatabase(lightECModel, 'display', [], humanAdapter);

%Show the differences between the model ec codes and the ones derived from the databases
sel = ~(cellfun(@isempty, lightECModel.ec.eccodes) & cellfun(@isempty, lightECModel2.ec.eccodes));
table(lightECModel.ec.eccodes(sel), lightECModel2.ec.eccodes(sel))

tic
%complexInfo = getComplexData('Homo sapiens', humanAdapter);
[lightECModel, foundComplex, proposedComplex] = applyComplexData(lightECModel, [], humanAdapter);
toc %Elapsed time is 62.289315 seconds - a bit slow, should be easy to optimize at some point.
%look at complex II - why does it not match in complexPortal?
%xxx = strcmp(lightECModel.ec.rxns, '001_MAR04652');
%protSel = lightECModel.ec.rxnEnzMat(xxx,:) > 0;
%sum(protSel) %4 - ok
%lightECModel.ec.enzymes(protSel)

kcatList = fuzzyKcatMatching(lightECModel, [], humanAdapter);
lightECModel = selectKcatValue(lightECModel,kcatList);
lightECModel = applyKcatConstraints(lightECModel);
constructEquations(lightECModel)
%hmm, look at the prot_pool usage
find(lightECModel.S(length(lightECModel.mets),:) < -10^12)
constructEquations(lightECModel, lightECModel.rxns(851))

lightECModel.ec.mw(strcmp(lightECModel.ec.rxns, '001_MAR03875')) %25195, reasonable
lightECModel.ec.kcat(strcmp(lightECModel.ec.rxns, '001_MAR03875')) %8.8300e-05, very small

%set protein pool constraint
lightECModel.ub(length(lightECModel.ub)) = 2.238315e-02 * 3600; %value taken from GECKO Light

lightECModelTuned = sensitivityTuning(lightECModel, 0.07, humanAdapter);




