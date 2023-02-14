%run this test case with the command
%results = runtests('geckoCoreFunctionTests.m')

function tests = geckoCoreFunctionTests
    tests = functiontests(localfunctions);
end

%code for writing the PhylDistStruct file for testing
%phylDistStruct = struct();
%phylDistStruct.ids = {'tst';'fls'};
%phylDistStruct.names = {'testus testus';'testus falsus'};
%phylDistStruct.distMat = [0 1; 1 0];
%save(fullfile(geckoPath,'test','unit_tests','ecTestGEM','data','PhylDist.mat'), 'phylDistStruct')


function testmakeEcModelFullModel_tc0001(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    
    %check that the ec model is
    %1) Expanded as expected
    %The reactions are sorted, so they come in alphabetical order here
    %Check rxns
    [expRxns,expMetNames,expS] = getBaseExpTstEcModelProperties();
    verifyEqual(testCase,ecModel.rxns,expRxns)
    %check mets
    verifyEqual(testCase,ecModel.metNames,expMetNames)
    %check S matrix
    verifyEqual(testCase,ecModel.S,expS)
    
    %2) Check the ec structure
    expEcRxns = {'R2_EXP_1';'R2_EXP_2';'R2_REV_EXP_1';'R2_REV_EXP_2';'R3';'R5'};
    verifyEqual(testCase,ecModel.ec.rxns,expEcRxns)
    verifyEqual(testCase,length(ecModel.ec.rxns),length(ecModel.ec.kcat)) % The kcats etc. are not set here, will be set later in the workflow
    verifyEqual(testCase,length(ecModel.ec.rxns),length(ecModel.ec.source))
    verifyEqual(testCase,length(ecModel.ec.rxns),length(ecModel.ec.notes))
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'}; %these are not set here, so should be empty - only check length
    verifyEqual(testCase,length(ecModel.ec.eccodes),length(expEccodes))
    expEcGenes = {'G1';'G2';'G3';'G4';'G5'};
    verifyEqual(testCase,ecModel.ec.genes,expEcGenes)
    expEnzymes = {'P1';'P2';'P3';'P4';'P5'};
    verifyEqual(testCase,ecModel.ec.enzymes,expEnzymes)
    expMW = [10000;20000;30000;40000;50000];
    verifyEqual(testCase,ecModel.ec.mw,expMW)
    expSequence = {'MRAL';'MNTD';'MSYN';'MDFM';'MLFK'};
    verifyEqual(testCase,ecModel.ec.sequence,expSequence);
    verifyEqual(testCase,length(expEcGenes),length(ecModel.ec.concs))
    expRxnEnzMat = sparse(length(expEcRxns),length(expEcGenes));
    expRxnEnzMat(1,1:2) = 1;
    expRxnEnzMat(2,3) = 1;
    expRxnEnzMat(3,1:2) = 1;
    expRxnEnzMat(4,3) = 1;
    expRxnEnzMat(5,4) = 1;
    expRxnEnzMat(6,5) = 1;
    verifyEqual(testCase,ecModel.ec.rxnEnzMat,full(expRxnEnzMat)) %TODO: should it be sparse of full?
end

function testmakeEcModelLightModel_tc0002(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, true, adapter);
    
    %check that the ec model is
    %1) Expanded as expected
    %The reactions are sorted, so they come in alphabetical order here
    %Check rxns
    expRxns = [model.rxns;'R1_REV';'R2_REV';'prot_pool_exchange'];
    verifyEqual(testCase,ecModel.rxns,expRxns)
    %check mets
    expMetNames = [model.metNames;'prot_pool'];
    verifyEqual(testCase,ecModel.metNames,expMetNames)
    %check S matrix
    expS = [model.S model.S(:,2:3)*-1 sparse(length(model.mets),1);sparse(1,length(ecModel.rxns)-1) 1];
    verifyEqual(testCase,ecModel.S,expS)
    
    %2) Check the ec structure
    expEcRxns = {'001_R2';'002_R2';'001_R3';'001_R5';'001_R2_REV';'002_R2_REV'};
    verifyEqual(testCase,ecModel.ec.rxns,expEcRxns)
    verifyEqual(testCase,length(ecModel.ec.rxns),length(ecModel.ec.kcat)) % The kcats etc. are not set here, will be set later in the workflow
    verifyEqual(testCase,length(ecModel.ec.rxns),length(ecModel.ec.source))
    verifyEqual(testCase,length(ecModel.ec.rxns),length(ecModel.ec.notes))
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3';'1.1.1.1';'1.1.1.1'}; %these are not set here, so should be empty - only check length
    verifyEqual(testCase,length(ecModel.ec.eccodes),length(expEccodes))
    expEcGenes = {'G1';'G2';'G3';'G4';'G5'};
    verifyEqual(testCase,ecModel.ec.genes,expEcGenes)
    expEnzymes = {'P1';'P2';'P3';'P4';'P5'};
    verifyEqual(testCase,ecModel.ec.enzymes,expEnzymes)
    expMW = [10000;20000;30000;40000;50000];
    verifyEqual(testCase,ecModel.ec.mw,expMW)
    expSequence = {'MRAL';'MNTD';'MSYN';'MDFM';'MLFK'};
    verifyEqual(testCase,ecModel.ec.sequence,expSequence);
    verifyEqual(testCase,length(expEcGenes),length(ecModel.ec.concs))
    expRxnEnzMat = sparse(length(expEcRxns),length(expEcGenes));
    expRxnEnzMat(1,1:2) = 1;
    expRxnEnzMat(2,3) = 1;
    expRxnEnzMat(3,4) = 1;
    expRxnEnzMat(4,5) = 1;
    expRxnEnzMat(5,1:2) = 1;
    expRxnEnzMat(6,3) = 1;
    verifyEqual(testCase,ecModel.ec.rxnEnzMat,full(expRxnEnzMat)) %TODO: should it be sparse of full?
end

function testapplyComplexDataFullModel_tc0003(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = applyComplexData(ecModel, [], adapter);
    
    expRxnEnzMat = sparse(6, 5);
    expRxnEnzMat(1,1:2) = [1 2];
    expRxnEnzMat(2,3) = 1;
    expRxnEnzMat(3,1:2) = [1 2];
    expRxnEnzMat(4,3) = 1;
    expRxnEnzMat(5,4) = 1;
    expRxnEnzMat(6,5) = 1;
    verifyEqual(testCase,ecModel.ec.rxnEnzMat,full(expRxnEnzMat)) %TODO: should it be sparse of full?
end

function testapplyComplexDataLightModel_tc0004(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, true, adapter);
    ecModel = applyComplexData(ecModel, [], adapter);
    
    expRxnEnzMat = sparse(6, 5);
    expRxnEnzMat(1,1:2) = [1 2];
    expRxnEnzMat(2,3) = 1;
    expRxnEnzMat(3,4) = 1;
    expRxnEnzMat(4,5) = 1;
    expRxnEnzMat(5,1:2) = [1 2];
    expRxnEnzMat(6,3) = 1;
    verifyEqual(testCase,ecModel.ec.rxnEnzMat,full(expRxnEnzMat)) %TODO: should it be sparse of full?
end

%For both full and light
function testsetProtPoolSize_tc0005(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = setProtPoolSize(ecModel, [], adapter);
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),10000)
    ecModel = setProtPoolSize(ecModel, 3);
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),3)

    %light
    ecModel = makeEcModel(model, true, adapter);
    ecModel = setProtPoolSize(ecModel, [], adapter);
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),10000)
    ecModel = setProtPoolSize(ecModel, 3);
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),3)
end

%For both full and light
function testgetECfromGEM_tc0006(testCase)
    %full
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel2 = getECfromGEM(ecModel);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromGEM(ecModel,ismember(ecModel.ec.rxns, {'R2_EXP_1';'R3'}));
    expEccodes = {'1.1.1.1';[];[];[];'1.1.2.1';[]};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    
    %light
    ecModel = makeEcModel(model, true, adapter);
    ecModel2 = getECfromGEM(ecModel);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3';'1.1.1.1';'1.1.1.1'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromGEM(ecModel,ismember(ecModel.ec.rxns, {'001_R2';'001_R3'}));
    expEccodes = {'1.1.1.1';[];'1.1.2.1';[];[];[]};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
end

%For both full and light
function testgetECfromDatabase_tc0007(testCase)
    %full
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel2 = getECfromDatabase(ecModel, 'display', [], adapter);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromDatabase(ecModel, 'display', ismember(ecModel.ec.rxns, {'R2_EXP_1';'R3'}), adapter);
    expEccodes = {'1.1.1.1';[];[];[];'1.1.2.1';[]};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    
    %light
    ecModel = makeEcModel(model, true, adapter);
    ecModel2 = getECfromDatabase(ecModel, 'display', [], adapter);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3';'1.1.1.1';'1.1.1.1'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromDatabase(ecModel, 'display', ismember(ecModel.ec.rxns, {'001_R2';'001_R3'}), adapter);
    expEccodes = {'1.1.1.1';[];'1.1.2.1';[];[];[]};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
end

function testModelAdapterManager_tc0008(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    verifyTrue(testCase,~isempty(adapter))
    ModelAdapterManager.setDefaultAdapter(adapter);
    verifyTrue(testCase,~isempty(ModelAdapterManager.getDefaultAdapter()))
    ModelAdapterManager.setDefaultAdapter([]);
    verifyTrue(testCase,isempty(ModelAdapterManager.getDefaultAdapter()))
    ModelAdapterManager.setDefaultAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    verifyTrue(testCase,~isempty(ModelAdapterManager.getDefaultAdapter()))
end

function testsaveECModel_tc0009(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModelFilledIn = saveECmodel(ecModel,'RAVEN','tmpTest','1',fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    loadedEcModel = readYAMLmodel(fullfile(geckoPath,'test','unit_tests','ecTestGEM','tmpTest.yml'), false);
    verifyEqual(testCase, ecModelFilledIn, loadedEcModel)
end


function testfuzzyKcatMatching_tc0010(testCase)
    %full
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapterFromPath(fullfile(geckoPath,'test','unit_tests','ecTestGEM'));
    model = getGeckoTestModel();
    %First all rxns, full model
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    kcatListLightFuzzy = fuzzyKcatMatching(ecModel, [], adapter);
    expKcats = [1;1;10;10;100;10];%substrate is more important than organism, and wildcard comes last
    verifyEqual(testCase,kcatListLightFuzzy.kcats, expKcats)
    expEcRxns = {'R2_EXP_1';'R2_EXP_2';'R2_REV_EXP_1';'R2_REV_EXP_2';'R3';'R5'};
    verifyEqual(testCase,kcatListLightFuzzy.rxns, expEcRxns)
    verifyEqual(testCase,kcatListLightFuzzy.substrates, {{'m1'};{'m1'};{'m2'};{'m2'};{'m1'};{'m2'}})
    verifyEqual(testCase,kcatListLightFuzzy.eccodes, {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'})
    verifyEqual(testCase,kcatListLightFuzzy.wildcardLvl, [0;0;0;0;1;1])
    verifyEqual(testCase,kcatListLightFuzzy.origin, [1;1;2;2;1;2])
    %Some rxns, full model
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    kcatListLightFuzzy = fuzzyKcatMatching(ecModel, ismember(ecModel.ec.rxns,{'R2_REV_EXP_1'}), adapter);
    expKcats = [10];%substrate is more important than organism, and wildcard comes last
    verifyEqual(testCase,kcatListLightFuzzy.kcats, expKcats)
    expEcRxns = {'R2_REV_EXP_1'};
    verifyEqual(testCase,kcatListLightFuzzy.rxns, expEcRxns)
    verifyEqual(testCase,kcatListLightFuzzy.substrates, {{'m2'}})
    verifyEqual(testCase,kcatListLightFuzzy.eccodes, {'1.1.1.1'})
    verifyEqual(testCase,kcatListLightFuzzy.wildcardLvl, [0])
    verifyEqual(testCase,kcatListLightFuzzy.origin, [2])
    
    %all rxns, light model
    ecModel = makeEcModel(model, true, adapter);
    ecModel = getECfromGEM(ecModel);
    kcatListLightFuzzy = fuzzyKcatMatching(ecModel, [], adapter);
    expKcats = [1;1;100;10;10;10];%substrate is more important than organism, and wildcard comes last
    verifyEqual(testCase,kcatListLightFuzzy.kcats, expKcats)
    expEcRxns = {'001_R2';'002_R2';'001_R3';'001_R5';'001_R2_REV';'002_R2_REV'};
    verifyEqual(testCase,kcatListLightFuzzy.rxns, expEcRxns)
    verifyEqual(testCase,kcatListLightFuzzy.substrates, {{'m1'};{'m1'};{'m1'};{'m2'};{'m2'};{'m2'}})
    verifyEqual(testCase,kcatListLightFuzzy.eccodes, {'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3';'1.1.1.1';'1.1.1.1'})
    verifyEqual(testCase,kcatListLightFuzzy.wildcardLvl, [0;0;1;1;0;0])
    verifyEqual(testCase,kcatListLightFuzzy.origin, [1;1;1;2;2;2])
    
    %some rxns, light model
    ecModel = makeEcModel(model, true, adapter);
    ecModel = getECfromGEM(ecModel);
    kcatListLightFuzzy = fuzzyKcatMatching(ecModel, ismember(ecModel.ec.rxns,{'001_R2_REV'}), adapter);
    expKcats = [10];%substrate is more important than organism, and wildcard comes last
    verifyEqual(testCase,kcatListLightFuzzy.kcats, expKcats)
    expEcRxns = {'001_R2_REV'};
    verifyEqual(testCase,kcatListLightFuzzy.rxns, expEcRxns)
    verifyEqual(testCase,kcatListLightFuzzy.substrates, {{'m2'}})
    verifyEqual(testCase,kcatListLightFuzzy.eccodes, {'1.1.1.1'})
    verifyEqual(testCase,kcatListLightFuzzy.wildcardLvl, [0])
    verifyEqual(testCase,kcatListLightFuzzy.origin, [2])
    
end




