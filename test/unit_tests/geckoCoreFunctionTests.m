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
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
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
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
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
    expS = [model.S model.S(:,2:3)*-1 sparse(length(model.mets),1);sparse(1,length(ecModel.rxns)-1) -1];
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
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = applyComplexData(ecModel, [], adapter, false);
    
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
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, true, adapter);
    ecModel = applyComplexData(ecModel, [], adapter, false);
    
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
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    verifyEqual(testCase,ecModel.lb(length(ecModel.rxns)),-1000)
    ecModel = setProtPoolSize(ecModel, 1, 5, 1);
    verifyEqual(testCase,ecModel.lb(length(ecModel.rxns)),-5000)

    %light
    ecModel = makeEcModel(model, true, adapter);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    verifyEqual(testCase,ecModel.lb(length(ecModel.rxns)),-1000)
    ecModel = setProtPoolSize(ecModel, 1, 5, 1);
    verifyEqual(testCase,ecModel.lb(length(ecModel.rxns)),-5000)
end

%For both full and light
function testgetECfromGEM_tc0006(testCase)
    %full
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel2 = getECfromGEM(ecModel);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromGEM(ecModel,ismember(ecModel.ec.rxns, {'R2_EXP_1';'R3'}));
    expEccodes = {'1.1.1.1';'';'';'';'1.1.2.1';''};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    
    %light
    ecModel = makeEcModel(model, true, adapter);
    ecModel2 = getECfromGEM(ecModel);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3';'1.1.1.1';'1.1.1.1'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromGEM(ecModel,ismember(ecModel.ec.rxns, {'001_R2';'001_R3'}));
    expEccodes = {'1.1.1.1';'';'1.1.2.1';'';'';''};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
end

%For both full and light
function testgetECfromDatabase_tc0007(testCase)
    %full
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel2 = getECfromDatabase(ecModel, [], 'display', adapter);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromDatabase(ecModel, ismember(ecModel.ec.rxns, {'R2_EXP_1';'R3'}),'display', adapter);
    expEccodes = {'1.1.1.1';'';'';'';'1.1.2.1';''};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    
    %light
    ecModel = makeEcModel(model, true, adapter);
    ecModel2 = getECfromDatabase(ecModel, [], 'display', adapter);
    expEccodes = {'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3';'1.1.1.1';'1.1.1.1'};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
    %for selected rxns only
    ecModel2 = getECfromDatabase(ecModel, ismember(ecModel.ec.rxns, {'001_R2';'001_R3'}), 'display', adapter);
    expEccodes = {'1.1.1.1';'';'1.1.2.1';'';'';''};
    verifyEqual(testCase,ecModel2.ec.eccodes,expEccodes)
end

function testModelAdapterManager_tc0008(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    verifyTrue(testCase,~isempty(adapter))
    ModelAdapterManager.setDefault(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    verifyTrue(testCase,~isempty(ModelAdapterManager.getDefault()))
    ModelAdapterManager.setDefault([]);
    verifyTrue(testCase,isempty(ModelAdapterManager.getDefault()))
end

function testsaveECModel_tc0009(testCase)
    % Test a round of model saving and loading
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    saveEcModel(ecModel, [], adapter);
    loadedEcModel = loadEcModel([],adapter);
    delete(fullfile(adapter.params.path,'models','ecModel.yml'));
    verifyEqual(testCase, ecModel, loadedEcModel)

    % Test loading of conventional GEM
    evalc('loadedModel = loadConventionalGEM('''',adapter);'); % Avoid throwing or warnings
    model=rmfield(model,{'annotation','date','description','version'});
    model.geneShortNames=model.genes;
    verifyEqual(testCase, model, loadedModel)
end

function testfuzzyKcatMatching_tc0010(testCase)
    %full
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    %First all rxns, full model
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    kcatListLightFuzzy = fuzzyKcatMatching(ecModel, [], adapter);
    expKcats = [1;1;10;10;100;1];%substrate is more important than organism, and wildcard comes last
    verifyEqual(testCase,kcatListLightFuzzy.kcats, expKcats)
    expEcRxns = {'R2_EXP_1';'R2_EXP_2';'R2_REV_EXP_1';'R2_REV_EXP_2';'R3';'R5'};
    verifyEqual(testCase,kcatListLightFuzzy.rxns, expEcRxns)
    verifyEqual(testCase,kcatListLightFuzzy.substrates, {{'m1'};{'m1'};{'m2'};{'m2'};{'m1'};{'m2'}})
    verifyEqual(testCase,kcatListLightFuzzy.eccodes, {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3'})
    verifyEqual(testCase,kcatListLightFuzzy.wildcardLvl, [0;0;0;0;1;1])
    verifyEqual(testCase,kcatListLightFuzzy.origin, [1;1;2;2;3;3])
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
    expKcats = [1;1;100;1;10;10];%substrate is more important than organism, and wildcard comes last
    verifyEqual(testCase,kcatListLightFuzzy.kcats, expKcats)
    expEcRxns = {'001_R2';'002_R2';'001_R3';'001_R5';'001_R2_REV';'002_R2_REV'};
    verifyEqual(testCase,kcatListLightFuzzy.rxns, expEcRxns)
    verifyEqual(testCase,kcatListLightFuzzy.substrates, {{'m1'};{'m1'};{'m1'};{'m2'};{'m2'};{'m2'}})
    verifyEqual(testCase,kcatListLightFuzzy.eccodes, {'1.1.1.1';'1.1.1.1';'1.1.2.1';'1.1.1.3';'1.1.1.1';'1.1.1.1'})
    verifyEqual(testCase,kcatListLightFuzzy.wildcardLvl, [0;0;1;1;0;0])
    verifyEqual(testCase,kcatListLightFuzzy.origin, [1;1;3;3;2;2])
    
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

%Tests mergeDLKcatAndFuzzyKcats, selectKcatValue, and applyKcatConstraints.
%Also to a certain extent tests writing of DLKcat files, but not the DLKCat algorithm or reading of the output
%In addition it tests that the small test model has the same growth rate for both full and light
function testKcats_tc0011(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    %add an extra R3 reaction to be able to check that wildcards go in if there is no kcat in the dlkcat list
    rxnsToAdd = struct();
    rxnsToAdd.rxns = {'R2a';'R3b'};
    rxnsToAdd.grRules = {'G1 and G2 or G3';'G4'};
    rxnsToAdd.equations = {'m1[c] => m2[c]'; 'm1[c] => m2[c]'};
    model = addRxns(model,rxnsToAdd, 3);
    model.eccodes{9} = '1.1.2.1';%no eccode for R2a, let dlkcat populate that
    
    %we only test with full, the model is not really involved in this code
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel = applyComplexData(ecModel, [], adapter, false);

    kcatListFuzzy = fuzzyKcatMatching(ecModel, [], adapter);
    %test to write a DLKcat
    filepath = fullfile(adapter.getParameters().path,'data','DLKcat_input_test.tsv');
    if exist(filepath, 'file')==2
      delete(filepath);
    end
    verifyTrue(testCase,~(exist(filepath, 'file')==2))
    [~, writtenTable] = evalc('writeDLKcatInput(ecModel, [], adapter, false, filepath)');
    verifyTrue(testCase,exist(filepath, 'file')==2) %check that we write a file, we don't check the contents
    if exist(filepath, 'file')==2 %clean up
      delete(filepath);
    end
    verifyEqual(testCase,writtenTable(1,:), {'R2_EXP_1','R2_EXP_1','R2_EXP_2','R2_REV_EXP_1','R2_REV_EXP_1','R2_REV_EXP_2','R2a_EXP_1','R2a_EXP_1','R2a_EXP_2','R3','R3b','R5'})
    verifyEqual(testCase,writtenTable(2,:), {'G1','G2','G3','G1','G2','G3','G1','G2','G3','G4','G4','G5'})
    verifyEqual(testCase,writtenTable(3,:), {'m1','m1','m1','m2','m2','m2','m1','m1','m1','m1','m1','m2'})
    %skip line 4, not set here since our fake metabolites don't have any smiles - could perhaps be fixed at some point
    verifyEqual(testCase,writtenTable(5,:), {'MRAL','MNTD','MSYN','MRAL','MNTD','MSYN','MRAL','MNTD','MSYN','MDFM','MDFM','MLFK'})
    %skip line 6, always set to 'NA' it seems
    
    %Create a suitable kcatlist from dlkcat
    dlkcatList = struct();
    dlkcatList.source = 'DLKcat';
    dlkcatList.rxns = {'R2_EXP_1';'R2_EXP_1';'R2_EXP_2';'R2_REV_EXP_1';'R2_REV_EXP_1';'R2_REV_EXP_2';'R2a_EXP_1';'R2a_EXP_1';'R2a_EXP_2';'R3';'R5'};
    dlkcatList.genes = {'G1';'G2';'G3';'G1';'G2';'G3';'G1';'G2';'G3';'G4';'G5'};
    dlkcatList.substrates = {'m1';'m1';'m1';'m2';'m2';'m2';'m1';'m1';'m1';'m1';'m2'};
    dlkcatList.kcats = [1001;1002;1003;1004;1005;1006;1007;1008;1009;1010;1011];
    mergedList = mergeDLKcatAndFuzzyKcats(dlkcatList, kcatListFuzzy, 6, 6, 1);%allow for use of wildcards
    
    %What we expect is that all R2, which have a good match (some with bad substrate) will be taken from fuzzy.
    %Furthermore, R3b will be taken from fuzzy, since we don't manage to predict it in dlkcat (we didn't include it in the dlkcat kcat list)
    %All R2a are missing the ec code, and therefore will be taken from dlkcat. 
    %R3 is a wildcard match with a value in the dlkcat list, and will thus be taken from dlkcat since it has higher prio
    %R5 doesn't have a kcat match in "brenda", and therefore uses dlkcat
    verifyEqual(testCase,mergedList.kcatSource, {'brenda';'brenda';'brenda';'brenda';'brenda';'DLKcat';'DLKcat';'DLKcat';'DLKcat';'DLKcat'})
    verifyEqual(testCase,mergedList.rxns, {'R2_EXP_1';'R2_EXP_2';'R2_REV_EXP_1';'R2_REV_EXP_2';'R3b';'R2a_EXP_1';'R2a_EXP_1';'R2a_EXP_2';'R3';'R5'})
    verifyEqual(testCase,mergedList.genes, {[];[];[];[];[];'G1';'G2';'G3';'G4';'G5'})
    verifyEqual(testCase,mergedList.substrates, {{'m1'};{'m1'};{'m2'};{'m2'};{'m1'};'m1';'m1';'m1';'m1';'m2'})
    verifyEqual(testCase,mergedList.kcats, [1;1;10;10;100;1007;1008;1009;1010;1011])
    verifyEqual(testCase,mergedList.eccodes, {'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.1.1';'1.1.2.1';[];[];[];[];[]})
    verifyEqual(testCase,mergedList.wildcardLvl, [0;0;0;0;1;NaN;NaN;NaN;NaN;NaN])
    verifyEqual(testCase,mergedList.origin, [1;1;2;2;3;NaN;NaN;NaN;NaN;NaN])%origin is 2 for testus falsus, 1 for the wildcard match which matches well on species and substrate

    %now test select
    %we expect the highest kcat value to be chosen in the R2a_EXP_1 case, i.e., use 1008, discard 1007
    ecModel = selectKcatValue(ecModel, mergedList);
    %{'R2_EXP_1';'R2_EXP_2';'R2_REV_EXP_1';'R2_REV_EXP_1';'R2_REV_EXP_2';'R2a_EXP_1';'R2a_EXP_2';'R3';'R3b';'R5'};
    expectedKcats = [1;1;10;10;1008;1009;1010;100;1011];
    verifyEqual(testCase,ecModel.ec.kcat, expectedKcats)
    
    %and apply - first full
    %Test a subset first
    ecModel = applyKcatConstraints(ecModel,{'R3';'R5'});
    %This should lead to a protein cost on reaction R3 (P4) and R5 (P5)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R1'))), [0;0;0;0;0],"AbsTol",10^-10) 
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R1_REV'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_EXP_1'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_EXP_2'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_REV_EXP_1'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_REV_EXP_2'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2a_EXP_1'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2a_EXP_2'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R3'))), [0;0;0;-40000/1010/3600;0],"AbsTol",10^-10) %MW 40000 (P4/G4), kcat 1010
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R3b'))), [0;0;0;0;0],"AbsTol",10^-10) %MW 40000 (P4/G4), kcat 1010
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R4'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R5'))), [0;0;0;0;-50000/1011/3600],"AbsTol",10^-10) %MW 50000 (P5/G5), kcat 1011
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'S1'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'S2'))), [0;0;0;0;0],"AbsTol",10^-10)
    %Now all
    ecModel = applyKcatConstraints(ecModel);
    
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R1'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R1_REV'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_EXP_1'))), [-10000/1/3600;-(2*20000)/1/3600;0;0;0],"AbsTol",10^-10) %MW 10000 + 2*20000 (P1 + 2*P2), kcat 1
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_EXP_2'))), [0;0;-30000/1/3600;0;0],"AbsTol",10^-10) %MW 30000 (P3), kcat 1
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_REV_EXP_1'))), [-10000/10/3600;-(2*20000)/10/3600;0;0;0],"AbsTol",10^-10) %MW 10000 + 2*20000 (P1 + 2*P2), kcat 10
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2_REV_EXP_2'))), [0;0;-30000/10/3600;0;0],"AbsTol",10^-10) %MW 30000 (P3), kcat 10
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2a_EXP_1'))), [-10000/1008/3600;-(2*20000)/1008/3600;0;0;0],"AbsTol",10^-10) %no ec code, so kcat from dlkcat, use the max kcat for both prot, i.e., max of 1007 and 1008
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R2a_EXP_2'))), [0;0;-30000/1009/3600;0;0],"AbsTol",10^-10) %no ec code, so kcat from dlkcat
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R3'))), [0;0;0;-40000/1010/3600;0],"AbsTol",10^-10) %MW 40000 (P4/G4), kcat 1010
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R3b'))), [0;0;0;-40000/100/3600;0],"AbsTol",10^-10) %MW 40000 (P4/G4), kcat 1010
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R4'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'R5'))), [0;0;0;0;-50000/1011/3600],"AbsTol",10^-10) %MW 50000 (P5/G5), kcat 1011
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'S1'))), [0;0;0;0;0],"AbsTol",10^-10)
    verifyEqual(testCase,full(ecModel.S(ismember(ecModel.mets, {'prot_P1';'prot_P2';'prot_P3';'prot_P4';'prot_P5'}),strcmp(ecModel.rxns, 'S2'))), [0;0;0;0;0],"AbsTol",10^-10)

    %Check getKcatAcrossIsozymes
    ecModel.ec.kcat(2)=0;
    ecModel=getKcatAcrossIsozymes(ecModel);
    verifyEqual(testCase,ecModel.ec.kcat, [1;1;10;10;1008;1009;1010;100;1011])

    %Check applyCustomKcats
    test = applyCustomKcats(ecModel, [], adapter);
    verifyEqual(testCase,test.ec.kcat, [100;200;50;50;100;200;1010;100;50])

    customKcats.proteins = {'P3'; 'P1 + P2'; ''};
    customKcats.kcat     = [200; 100; 50];
    customKcats.rxns     = {'';'';'R2_REV, R5'};

    test = applyCustomKcats(ecModel, customKcats, adapter);
    verifyEqual(testCase,test.ec.kcat, [100;200;50;50;100;200;1010;100;50])


    %now apply for light
    %%%%%%%%%%%%%%%%%%%
    
    lecModel = makeEcModel(model, true, adapter);
    lecModel = getECfromGEM(lecModel);
    lecModel = applyComplexData(lecModel, [], adapter, false);
    kcatListFuzzy = fuzzyKcatMatching(lecModel, [], adapter);
    
    %Create a suitable kcatlist from dlkcat
    dlkcatList = struct();
    dlkcatList.source = 'DLKcat';
    dlkcatList.rxns = {'001_R2';'001_R2';'002_R2';'001_R3';'001_R5';'001_R2a';'001_R2a';'002_R2a';'001_R2_REV';'001_R2_REV';'002_R2_REV'};
    dlkcatList.genes = {'G1';'G2';'G3';'G4';'G5';'G1';'G2';'G3';'G1';'G2';'G3'};
    dlkcatList.substrates = {'m1';'m1';'m1';'m1';'m2';'m1';'m1';'m1';'m2';'m2';'m2'};
    dlkcatList.kcats = [1001;1002;1003;1010;1011;1007;1008;1009;1004;1005;1006];

    mergedList = mergeDLKcatAndFuzzyKcats(dlkcatList, kcatListFuzzy, 6, 6, 1);%allow for use of wildcards
    lecModel = selectKcatValue(lecModel, mergedList);
    
    %and apply - first full
    %Test a subset first
    lecModel = applyKcatConstraints(lecModel,{'R3';'R5'});
    %This should lead to a protein cost on reaction R3 (P4) and R5 (P5)
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R1'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2_REV'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2a'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R3'))),-40000/1010/3600,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R3b'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R4'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R5'))),-50000/1011/3600,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'S1'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'S2'))),0,"AbsTol",10^-10) 
    %Now all
    lecModel = applyKcatConstraints(lecModel);
    
    %This should be the same as for full, with the difference that the R2 reactions are just one reaction where the minimum cost is chosen among isozymes
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R1'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2'))),-min(10000/1/3600 + (2*20000)/1/3600, 30000/1/3600),"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2_REV'))),-min(10000/10/3600 + (2*20000)/10/3600, 30000/10/3600),"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2a'))),-min(10000/1008/3600+(2*20000)/1008/3600, 30000/1009/3600),"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R3'))),-40000/1010/3600,"AbsTol",10^-10)
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R3b'))),-40000/100/3600,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R4'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R5'))),-50000/1011/3600,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'S1'))),0,"AbsTol",10^-10) 
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'S2'))),0,"AbsTol",10^-10) 
    
    %Finally test if full and light has the same growth rate
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    lecModel = setProtPoolSize(lecModel, [], [], [], adapter);
    
    res = solveLP(ecModel,1);
    lres = solveLP(lecModel,1);
    verifyEqual(testCase,res.f,lres.f,"AbsTol",10^-10) 
end

%Does not test the data download, only operates from a stored file
function testfindMetSmiles_tc0012(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    [~, ecModel] = evalc("findMetSmiles(ecModel, adapter, false)");

    verifyEqual(testCase,ecModel.metSmiles,{'C(C1C)O';'C1C(=NC2)';'C(C1C)O';'C1C(=NC2)';'';'';'';'';'';''})
end

%Tests readProteomics, constrainEnzConcs, flexibilizeEnzConcs, and getConcControlCoeffs.
function testProteomcisIntegration_tc0013(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel = applyComplexData(ecModel, [], adapter, false);
    kcatListFuzzy = fuzzyKcatMatching(ecModel, [], adapter);
    ecModel  = selectKcatValue(ecModel, kcatListFuzzy);
    ecModel  = applyKcatConstraints(ecModel);
    ecModel  = setProtPoolSize(ecModel,[],[],[],adapter);

    % test that proteomics data is correct loaded into protData
    protData = loadProtData(1,[],[],adapter);
    verifyEqual(testCase,protData.abundances,[0.7292388;0.03692241;0.318175;5.1959184;0.15647268])
    verifyEqual(testCase,protData.uniprotIDs,{'P1';'P2';'P3';'P4';'P5'})

    % test to protData is correct included in the model
    ecModel = fillEnzConcs(ecModel,protData);
    verifyEqual(testCase,ecModel.ec.concs,[0.7292388;0.03692241;0.318175;5.1959184;0.15647268])

    % test that usage protein are correctly constraint
    [~, usageRxnIdx] = ismember(strcat('usage_prot_', ecModel.ec.enzymes), ecModel.rxns);
    ecModel = constrainEnzConcs(ecModel);
    verifyEqual(testCase,ecModel.lb(usageRxnIdx),-ecModel.ec.concs)

    % test that usage protein are correctly constraint. Sol.f give 0.1127,
    % increse objective up to 0.5
    [~, ecModel, flexEnz] =  evalc("flexibilizeEnzConcs(ecModel, 0.4,[],[],adapter,false)");
    [~, usageRxnIdx] = ismember(strcat('usage_prot_', flexEnz.uniprotIDs), ecModel.rxns);
    verifyEqual(testCase,ecModel.lb(usageRxnIdx),-flexEnz.flexConcs)
end

