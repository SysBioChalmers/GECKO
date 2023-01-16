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

kcatList = fuzzyKcatMatching(lightECModel, [], humanAdapter);

ligthECModelWithoutEC = lightECModel;
ligthECModelWithoutEC.ec.eccodes = [];
ltECWithNewEC = getECfromDatabase(ligthECModelWithoutEC, 'display', [], humanAdapter);%Doesn't work so well, conflicts in uniprot db




%just fill the light model with random kcats for now
lightECModelFilledIn = lightECModel;
lightECModelFilledIn.ec.kcat = unifrnd(2,10^5, length(lightECModel.ec.kcat),1);

lightECModelFilledIn = applyKcatConstraints(lightECModelFilledIn);
lightECModelFilledIn.S(strcmp(lightECModelFilledIn.mets, 'prot_pool'), :)

constructEquations(lightECModelFilledIn)





%{
x = strlength(modSimp.grRules);
x(x > 1000)

x2 = strlength(modExp.grRules);

table(x(x2> x), x2(x2>x))
modSimp.grRules(x2> x)
modExp.grRules(x2> x)


modSimp.rxns(x > 1000)
modSimp.grRules(x > 1000)
constructEquations(modSimp, modSimp.rxns(x > 1000))
modSimp.eccodes(x > 1000)

modExp.grRules(strcmp(modExp.rxns, 'MAR04137'))
%}


%{
tic
rev = convertToIrrev(modExp);
toc%Elapsed time is 0.105283 seconds.

tic
expanded = expandModel(rev);
toc%Elapsed time is 419.827173 seconds.

tic
expanded2 = expandModel2(rev);
toc %Elapsed time is 22.146821 seconds.

all(all(expanded2.S == expanded.S))%ok
all(all(expanded2.rxnGeneMat == expanded.rxnGeneMat))%ok
all(strcmp(expanded2.rxns,expanded.rxns))%ok
all(strcmp(expanded2.rxnNames,expanded.rxnNames))%ok
all(expanded2.ub == expanded.ub)%ok
all(expanded2.lb == expanded.lb)%ok
all(expanded2.rev == expanded.rev)%ok
all(expanded2.c == expanded.c)%ok
res = true;
for i = 1:length(expanded2.subSystems)
    res = res & all(strcmp(expanded2.subSystems{i},expanded.subSystems{i}));
end
res %ok
all(strcmp(expanded2.eccodes,expanded.eccodes))%ok
all(strcmp(expanded2.rxnFrom,expanded.rxnFrom))%ok
all(strcmp(expanded2.rxnNotes,expanded.rxnNotes))%ok
all(strcmp(expanded2.rxnReferences,expanded.rxnReferences))%ok
all(expanded2.rxnConfidenceScores == expanded.rxnConfidenceScores)%ok
all(strcmp(expanded2.grRules,expanded.grRules))%ok

S2 = [rev.S rev.S(:,randsample(size(rev.S,2),size(rev.S,2)))];

%only .rxns
%Elapsed time is 6.537047 seconds.

%before .S
%Elapsed time is 21.530236 seconds.

profile on  -history
expanded2 = expandModel2(rev);
p = profile('info')

numEvents = size(p.FunctionHistory,2);
for n = 1:numEvents
    name = p.FunctionTable(p.FunctionHistory(2,n)).FunctionName;
    
    if p.FunctionHistory(1,n) == 0
        disp(['Entered ' name]);
    else
        disp(['Exited ' name]);
    end
end

table({p.FunctionTable.FunctionName}.', {p.FunctionTable.TotalTime}.')
profile off

toc
%}


%We might want to add a call to standardizeGrRules in makeEcModel.
%We could then remove the warning code in expandModel.

