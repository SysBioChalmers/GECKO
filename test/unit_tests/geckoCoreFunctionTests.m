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
    %check S matrix (forward prot_pool_exchange: -> prot_pool, coeff +1)
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
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),1000)
    ecModel = setProtPoolSize(ecModel, 1, 5, 1);
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),5000)

    %light
    ecModel = makeEcModel(model, true, adapter);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),1000)
    ecModel = setProtPoolSize(ecModel, 1, 5, 1);
    verifyEqual(testCase,ecModel.ub(length(ecModel.rxns)),5000)
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
    verifyEqual(testCase,ecModel.ub(usageRxnIdx),ecModel.ec.concs)

    % test that usage protein are correctly constraint. Sol.f give 0.1127,
    % increse objective up to 0.5
    [~, ecModel, flexEnz] =  evalc("flexibilizeEnzConcs(ecModel, 0.4,[],[],adapter,false)");
    [~, usageRxnIdx] = ismember(strcat('usage_prot_', flexEnz.uniprotIDs), ecModel.rxns);
    verifyEqual(testCase,ecModel.ub(usageRxnIdx),flexEnz.flexConcs)
end

function testCalculateFfactor_tc0014(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);

    % ecTestGEM ships no data/paxDB.tsv, so with no protData supplied calculateFfactor must
    % fall back to its documented default of 0.5 --- previously this crashed instead,
    % since execution fell through past the default-0.5 branch and went on to index into
    % protData as though it were the struct that branch never produced.
    [~, f] = evalc("calculateFfactor(ecModel, [], [], adapter)");
    verifyEqual(testCase,f,0.5)
end

function testWriteDLKcatInputSubset_tc0015(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    rxnsToAdd = struct();
    rxnsToAdd.rxns = {'R2a'};
    rxnsToAdd.grRules = {'G1 and G2 or G3'};
    rxnsToAdd.equations = {'m1[c] <=> m2[c]'};
    model = addRxns(model,rxnsToAdd, 3); %no eccode for R2a, so it never gets a fuzzy match

    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);

    % A non-trivial ecRxns subset --- neither "all reactions" (ecRxns=[]) nor a subset
    % that happens to start at index 1 --- is the ordinary, documented way to call
    % writeDLKcatInput (e.g. "only the reactions fuzzy matching couldn't find a kcat
    % for"), but was never exercised by any existing test: every prior call site here
    % and in the manual tests either omits ecRxns or passes a full-length mask. With
    % R2a_EXP_1/R2a_EXP_2 in the middle of ec.rxns rather than at its start, this used
    % to silently write zero rows (an internal filtering step re-used already-consumed
    % indices from ec.rxns space to index into the already-filtered, subset-length
    % array, clearing every row it should have kept), and after a first, incomplete fix
    % it wrote the right *number* of rows but attributed them to the wrong reactions
    % (the same re-used-index mistake, one step further down).
    ecRxns = ismember(ecModel.ec.rxns, {'R2a_EXP_1','R2a_EXP_2'});
    filepath = fullfile(adapter.getParameters().path,'data','DLKcat_input_subset_test.tsv');
    if exist(filepath, 'file')==2
      delete(filepath);
    end
    writtenTable = writeDLKcatInput(ecModel, 'ecRxns', ecRxns, 'modelAdapter', adapter, ...
        'onlyWithSmiles', false, 'filename', filepath, 'overwrite', true);
    if exist(filepath, 'file')==2 %clean up
      delete(filepath);
    end
    verifyEqual(testCase,writtenTable(1,:), {'R2a_EXP_1','R2a_EXP_1','R2a_EXP_2'})
    verifyEqual(testCase,writtenTable(2,:), {'G1','G2','G3'})
    verifyEqual(testCase,writtenTable(3,:), {'m1','m1','m1'})
end


function testTestGEMAdapterSpontaneousReactions_tc0016(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));

    % getSpontaneousReactions used to reference an undefined variable
    % (rxns_tsv.rxns) and crash unconditionally whenever called -- this is the
    % only method that ever calls it (from getStandardKcat.m), so
    % getStandardKcat could never run against TestGEMAdapter at all.
    model = getGeckoTestModel();
    [spont, spontRxnNames] = adapter.getSpontaneousReactions(model);
    verifyEqual(testCase,find(spont),5)
    verifyEqual(testCase,spontRxnNames,{'R4'})

    % Fixing the undefined variable alone was not enough: the position it set
    % (spont(5) = true) is only valid for this 7-reaction conventional model.
    % getStandardKcat's only real caller passes the already-expanded ecModel,
    % where R4's position among model.rxns has moved -- matching by reaction
    % id instead of position (mirroring geckopy's own TestGEMAdapter port)
    % survives that.
    ecModel = makeEcModel(model, false, adapter);
    [spontEc, spontRxnNamesEc] = adapter.getSpontaneousReactions(ecModel);
    verifyEqual(testCase,find(spontEc),find(strcmp(ecModel.rxns,'R4')))
    verifyEqual(testCase,spontRxnNamesEc,{'R4'})
end


function testApplyKcatConstraintsDuplicateAccession_tc0017(testCase)
    % Two ec.enzymes rows can map to the same prot_<accession> metabolite --
    % distinct genes sharing one UniProt entry. applyKcatConstraints used to
    % write model.S(linearIndices) = -newKcats(:,4) directly: a plain indexed
    % assignment on a repeated linear index keeps only the last value written,
    % silently dropping one enzyme's contribution instead of adding it to the
    % other's.
    clear model
    model.rxns = {'R1'};
    model.mets = {'prot_P1'};
    model.S = sparse(1,1);
    model.lb = -1000; model.ub = 1000;
    model.ec.geckoLight = false;
    model.ec.rxns = {'R1'};
    model.ec.kcat = 100;
    model.ec.enzymes = {'P1'; 'P1'};
    model.ec.mw = [10000; 20000];
    model.ec.rxnEnzMat = sparse([1 1]);

    model = applyKcatConstraints(model);
    expected = -(10000+20000)/100/3600; %summed, not just the second row's -20000/100/3600
    verifyEqual(testCase,full(model.S(1,1)),expected,'AbsTol',1e-12)
end


function testAddNewRxnsToEC_tc0018(testCase)
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);

    % Two new reactions, each with its own isozyme (OR) grRule, added in the same call.
    % Regression test: addNewRxnsToEC used to lose track of which reaction it was
    % splitting once more than one needed splitting. It removed (and appended) one entry
    % at a time while iterating a list of indices computed before any removal, so each
    % removal shifted the positions of every later entry --- the second (and any further)
    % reaction's grRule was then read from whatever had shifted into its old position,
    % rather than its own.
    newRxns.rxns      = {'RNEWA'; 'RNEWB'};
    newRxns.rxnNames  = {'RNEWA'; 'RNEWB'};
    newRxns.equations = {'m1[c] => e2[e]'; 'm2[c] => e1[e]'};
    newRxns.grRules   = {'GNEW1 or GNEW2'; 'GNEW3 or GNEW4'};

    newEnzymes.enzymes = {'ENEW1'; 'ENEW2'; 'ENEW3'; 'ENEW4'};
    newEnzymes.genes   = {'GNEW1'; 'GNEW2'; 'GNEW3'; 'GNEW4'};
    newEnzymes.mw      = [15000; 25000; 35000; 45000];

    [ecModel, rxnsAdded, enzAdded] = addNewRxnsToEC(ecModel, newRxns, newEnzymes, adapter);

    expRxnsAdded = {'RNEWA_EXP_1'; 'RNEWA_EXP_2'; 'RNEWB_EXP_1'; 'RNEWB_EXP_2'};
    verifyEqual(testCase, sort(rxnsAdded), sort(expRxnsAdded))
    verifyEqual(testCase, sort(enzAdded), sort(newEnzymes.enzymes))

    % Each split reaction must carry exactly the one gene it was split from --- not the
    % unsplit "X or Y" string the bug left on the second (and later) reaction, nor a
    % double-suffixed id from re-splitting an already-split entry.
    [~, idxA1] = ismember('RNEWA_EXP_1', ecModel.rxns);
    [~, idxA2] = ismember('RNEWA_EXP_2', ecModel.rxns);
    [~, idxB1] = ismember('RNEWB_EXP_1', ecModel.rxns);
    [~, idxB2] = ismember('RNEWB_EXP_2', ecModel.rxns);
    verifyEqual(testCase, ecModel.grRules{idxA1}, 'GNEW1')
    verifyEqual(testCase, ecModel.grRules{idxA2}, 'GNEW2')
    verifyEqual(testCase, ecModel.grRules{idxB1}, 'GNEW3')
    verifyEqual(testCase, ecModel.grRules{idxB2}, 'GNEW4')

    % And the enzyme-coupling matrix must reflect the same one-gene-per-split-reaction
    % mapping.
    [~, ecA1] = ismember('RNEWA_EXP_1', ecModel.ec.rxns);
    [~, ecA2] = ismember('RNEWA_EXP_2', ecModel.ec.rxns);
    [~, ecB1] = ismember('RNEWB_EXP_1', ecModel.ec.rxns);
    [~, ecB2] = ismember('RNEWB_EXP_2', ecModel.ec.rxns);
    [~, enz1] = ismember('ENEW1', ecModel.ec.enzymes);
    [~, enz2] = ismember('ENEW2', ecModel.ec.enzymes);
    [~, enz3] = ismember('ENEW3', ecModel.ec.enzymes);
    [~, enz4] = ismember('ENEW4', ecModel.ec.enzymes);
    verifyEqual(testCase, ecModel.ec.rxnEnzMat(ecA1, enz1), 1)
    verifyEqual(testCase, ecModel.ec.rxnEnzMat(ecA2, enz2), 1)
    verifyEqual(testCase, ecModel.ec.rxnEnzMat(ecB1, enz3), 1)
    verifyEqual(testCase, ecModel.ec.rxnEnzMat(ecB2, enz4), 1)
end


function testReportEnzymeUsageTopAbsUsageOutOfBounds_tc0019(testCase)
    % reportEnzymeUsage's default topAbsUsage (10) used to reach
    % usageData.protID(topUse(1:topAbsUsage)) unclamped: on any model with fewer
    % than 10 enzymes -- ecTestGEM has 5 -- that indexed past the end of topUse
    % and crashed. Regression test for the default call, plus the two other
    % "all enzymes" spellings the docstring promises (0 and Inf), plus an
    % ordinary in-range value left unaffected.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 100;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);

    fluxes = zeros(numel(ecModel.rxns),1);
    usageData = enzymeUsage(ecModel, fluxes);
    nEnz = numel(usageData.protID);

    report = reportEnzymeUsage(ecModel, usageData);
    verifyEqual(testCase, height(report.topAbsUsage), nEnz)

    report = reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', 0);
    verifyEqual(testCase, height(report.topAbsUsage), nEnz)

    report = reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', Inf);
    verifyEqual(testCase, height(report.topAbsUsage), nEnz)

    report = reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', 2);
    verifyEqual(testCase, height(report.topAbsUsage), 2)
end


function testSigmaFitterModelMatchesReturnedSigma_tc0020(testCase)
    % sigmaFitter's own docstring promises the returned model has its protein
    % pool "adapted to the optimal sigma-factor" -- but the grid-search loop
    % left `model` sized for whichever sigma the loop tried *last* (i=100,
    % i.e. sigma=1), not the best-fitting one `sigma` reports. Regression test:
    % the returned model's protein-pool upper bound must equal Ptot*f*sigma*1000
    % at the *returned* sigma, not at sigma=1 (unless they happen to coincide).
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns,'R3')) = 5;
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns,'R5')) = 200;
    ecModel.ec.source(:) = {'manual'};
    ecModel.c = double(strcmp(ecModel.rxns, 'R5'));
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);

    [fittedModel, sigma] = sigmaFitter(ecModel, 'growthRate', 100, 'Ptot', 0.5, ...
        'f', 4, 'makePlot', false, 'modelAdapter', adapter);

    % This growth target does not resolve to sigma=1 on this fixture, so a model
    % left at the loop's last trial is a real, detectable mismatch, not a
    % coincidental pass.
    verifyNotEqual(testCase, sigma, 1)

    poolIdx = strcmp(fittedModel.rxns, 'prot_pool_exchange');
    expected = 0.5 * 4 * sigma * 1000;
    verifyEqual(testCase, fittedModel.ub(poolIdx), expected, 'AbsTol', 1e-9)
end


function testFuzzyKcatMatchingWildcardIsPrefixNotSubstring_tc0021(testCase)
    % fuzzyKcatMatching's wildcard escalation truncates an EC code to a prefix (e.g.
    % '1.1.1.1' -> '1.-.-.-' -> truncated query '1.') and searches BRENDA's EC column
    % for it via extract_string_matches. The optimized path (used whenever the
    % ECIndexIds/EcIndexIndices index is available, i.e. always, from within
    % fuzzyKcatMatching itself) used contains(...) instead of an anchored prefix
    % match, so a query for enzyme class 1 ('1.') could also match an EC code from an
    % unrelated class that merely happens to contain the same two characters
    % somewhere else in its string -- '4.2.1.1' contains '1.' between its third and
    % fourth levels.
    %
    % The fixture has one legitimate class-1 entry (a real wildcard match, giving the
    % escalation loop something to succeed on and stop at) and one unrelated
    % '4.2.1.1' entry with a deliberately *larger* kcat: mainMatch takes the max
    % among same-level matches, so the bug (matching both) picks the larger, wrong
    % value; the fix (matching only the real class-1 entry) picks the smaller, correct
    % one.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);

    scratchDir = fullfile(tempname);
    brendaDir = fullfile(scratchDir, 'data');
    mkdir(brendaDir);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    % Organism matches TestGEMAdapter's own org_name exactly, so the "any organism"
    % match tier resolves trivially instead of falling through phylDist's
    % KEGG/genus lookup (which a fabricated organism name cannot resolve, causing
    % every candidate to be filtered out regardless of the EC match itself).
    fid = fopen(fullfile(brendaDir, 'max_KCAT.txt'), 'w');
    fprintf(fid, 'EC1.9.9.9\tm1\ttestus testus//*//*\t42\t*\n');
    fprintf(fid, 'EC4.2.1.1\tm1\ttestus testus//*//*\t500\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'max_MW.txt'), 'w');
    fprintf(fid, 'EC0.0.0.0\t*\tplaceholder//*//*\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'max_SA.txt'), 'w');
    fprintf(fid, 'EC0.0.0.0\t*\tplaceholder//*//*\t1\t*\n');
    fclose(fid);
    % getPhylDistStructPath also resolves under params.path/data; reuse ecTestGEM's
    % own real fixture rather than fabricating a second one.
    copyfile(fullfile(geckoPath,'test','unit_tests','ecTestGEM','data','PhylDist.mat'), ...
        fullfile(brendaDir, 'PhylDist.mat'));

    % TestGEMAdapter.getBrendaDBFolder resolves to <params.path>/data --- value class,
    % so redirecting a copy cannot affect the original (same pattern kcat_chain_ectestgem
    % and other scenarios use for scenario-specific data).
    fuzzyAdapter = adapter;
    fuzzyAdapter.params.path = scratchDir;

    kcatList = fuzzyKcatMatching(ecModel, ismember(ecModel.ec.rxns,{'R2_EXP_1'}), fuzzyAdapter);
    verifyEqual(testCase, kcatList.kcats, 42)
end


function testFuzzyKcatMatchingTieBreakRespectsWildcardCount_tc0022(testCase)
    % iterativeMatch's cross-EC-token tie-break is supposed to (1) keep only the EC
    % tokens that matched with the fewest wildcards, (2) among those keep the ones
    % with the best (lowest) origin, then (3) take the largest kcat. Step 2's
    % `best_pos = (origin == min(new_origin(new_origin~=0)))` re-evaluates
    % `origin == ...` against the *full*, unfiltered origin vector instead of the
    % wildcard-filtered subset from step 1 -- so any EC token that matched at a
    % *worse* (higher) wildcard level, but happens to share the same origin value,
    % is let back in. Step 3 then maximises kcat across that too-permissive set, so
    % a less-specific, more-wildcarded match can silently win over a fully-specific
    % one whenever it happens to report a bigger kcat.
    %
    % The fixture gives R2_EXP_1 two EC tokens: '9.9.9.1' matches BRENDA directly
    % (0 wildcards, kcat 42); '9.9.8.-' only matches after escalating twice more (2
    % wildcards, kcat 500) -- both land on origin 3 (organism matches, substrate
    % deliberately does not, so origin 1/2 are skipped). The correct answer is the
    % 0-wildcard match (42); the bug picks the 2-wildcard one (500) because both
    % share origin 3.
    %
    % Found while running the kcat-gathering pipeline against real BRENDA data and
    % yeast-GEM at genome scale in raven-gecko-parity: ~119 reactions got a kcat
    % that genuinely differed (not just missing) between the MATLAB and Python
    % implementations even after two other confirmed bugs were fixed; this
    % tie-break scoping bug, triggered whenever a multi-EC-code reaction's tokens
    % match at different wildcard levels but the same origin, accounts for at
    % least some of those.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.eccodes{strcmp(ecModel.ec.rxns, 'R2_EXP_1')} = '9.9.9.1;9.9.8.-';

    scratchDir = fullfile(tempname);
    brendaDir = fullfile(scratchDir, 'data');
    mkdir(brendaDir);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    % Substrate deliberately does not match the model's actual substrate (m1), so
    % origin 1/2 (which require a substrate match) are skipped and both tokens are
    % forced to resolve at origin 3, isolating the wildcard-count tie-break.
    fid = fopen(fullfile(brendaDir, 'max_KCAT.txt'), 'w');
    fprintf(fid, 'EC9.9.9.1\tunrelatedsubstrate\ttestus testus//*//*\t42\t*\n');
    fprintf(fid, 'EC9.9.7.5\tunrelatedsubstrate\ttestus testus//*//*\t500\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'max_MW.txt'), 'w');
    fprintf(fid, 'EC0.0.0.0\t*\tplaceholder//*//*\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'max_SA.txt'), 'w');
    fprintf(fid, 'EC0.0.0.0\t*\tplaceholder//*//*\t1\t*\n');
    fclose(fid);
    % getPhylDistStructPath also resolves under params.path/data; reuse ecTestGEM's
    % own real fixture rather than fabricating a second one.
    copyfile(fullfile(geckoPath,'test','unit_tests','ecTestGEM','data','PhylDist.mat'), ...
        fullfile(brendaDir, 'PhylDist.mat'));

    % TestGEMAdapter.getBrendaDBFolder resolves to <params.path>/data --- value class,
    % so redirecting a copy cannot affect the original (same pattern kcat_chain_ectestgem
    % and other scenarios use for scenario-specific data).
    fuzzyAdapter = adapter;
    fuzzyAdapter.params.path = scratchDir;

    kcatList = fuzzyKcatMatching(ecModel, ismember(ecModel.ec.rxns,{'R2_EXP_1'}), fuzzyAdapter);
    verifyEqual(testCase, kcatList.kcats, 42)
    verifyEqual(testCase, kcatList.wildcardLvl, 0)
end


function testWriteOpenKineticsPredictorInputReactionSubsetIndex_tc0023(testCase)
    % writeOpenKineticsPredictorInput restricts its work to the requested ecRxns
    % subset by first building clearedS = reducedS(:, origRxnIdxs), then finding
    % substrates with reactionIdxs = find(clearedS < 0) -- reactionIdxs is a
    % *column position within that subset* (1..numel(ecRxnsIdx)), not an absolute
    % index into model.ec.rxns. The next line indexed model.ec.rxnEnzMat directly
    % with reactionIdxs, silently picking the wrong ec.rxns row (or the wrong
    % protein set entirely) whenever the requested subset excludes any earlier
    % ec.rxns entry -- which shifts every later entry's subset-position below its
    % absolute position.
    %
    % ecTestGEM's R2 has two isozymes (G1+G2 complex, and G3 alone), each
    % expanded to its own ec.rxns row (R2_EXP_1, R2_EXP_2) plus a reverse copy
    % (R2_REV_EXP_1, R2_REV_EXP_2) -- ec.rxns order is [R2_EXP_1, R2_EXP_2,
    % R2_REV_EXP_1, R2_REV_EXP_2, R3, R5]. Requesting ecRxns={R2_EXP_1, R2_EXP_2,
    % R3} excludes the two REV entries in between, so R3's subset position (3)
    % no longer equals its absolute ec.rxns position (5): the bug reads
    % rxnEnzMat row 3 (R2_REV_EXP_1's genes, G1+G2) instead of row 5 (R3's own
    % gene, G4), so R3/G4 never appears in the output and G1/G2 appear twice.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);

    % Self-contained SMILES: only m1c (R2/R3's shared substrate) needs one for
    % this test, so skip findMetSmiles/PubChem entirely.
    ecModel.metSmiles = repmat({''}, numel(ecModel.mets), 1);
    ecModel.metSmiles(strcmp(ecModel.mets,'m1c')) = {'C(C1C)O'};

    scratchDir = fullfile(tempname);
    mkdir(fullfile(scratchDir,'data'));
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    % writeOpenKineticsPredictorInput resolves its own output path through
    % params.path/data --- value-class copy, same redirect pattern as every
    % other scenario/test in this file that writes a file.
    writeAdapter = adapter;
    writeAdapter.params.path = scratchDir;

    writeOpenKineticsPredictorInput(ecModel, 'ecRxns', ismember(ecModel.ec.rxns, {'R2_EXP_1','R2_EXP_2','R3'}), ...
        'modelAdapter', writeAdapter, 'onlyWithSmiles', true, 'overwrite', true);

    fID = fopen(fullfile(scratchDir,'data','OKP.csv'));
    raw = textscan(fID, '%s %s', 'Delimiter', ',', 'HeaderLines', 1);
    fclose(fID);
    rows = sort(strcat(raw{1}, ',', raw{2}));

    expected = sort({'MRAL,C(C1C)O'; 'MNTD,C(C1C)O'; 'MSYN,C(C1C)O'; 'MDFM,C(C1C)O'});
    verifyEqual(testCase, rows, expected)
end


function testGetEnzymeBottlenecksRanksByShadowPrice_tc0024(testCase)
    % getEnzymeBottlenecks: solves the model, ranks enzymes by |shadow price| of
    % their prot_<id> mass-balance constraint, returns the top N. Cross-verified
    % against geckopy's get_enzyme_bottlenecks on the identical fixture (uniform
    % kcat=10 across every ec.rxns row, prot pool sized via Ptot=0.5/f=0.5/
    % sigma=0.5): both sides agree exactly on the objective (90), every column
    % (shadowPrice=-0.72 for all five enzymes, flux/capUsage/upperBound), and the
    % top-3 selection (P1/P2/P3, tied at -0.72, in original ec.enzymes order).
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 10;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, 0.5, 0.5, 0.5, adapter);

    sol = solveLP(ecModel);
    verifyEqual(testCase, sol.stat, 1)
    verifyEqual(testCase, sol.f, 90, 'AbsTol', 1e-9)

    bottlenecks = getEnzymeBottlenecks(ecModel);
    verifyEqual(testCase, height(bottlenecks), 5)
    verifyEqual(testCase, bottlenecks.shadowPrice, repmat(-0.72, 5, 1), 'AbsTol', 1e-9)
    [~, order] = ismember({'P1','P2','P3','P4','P5'}, bottlenecks.uniprot);
    verifyEqual(testCase, bottlenecks.flux(order), [0;0;0;0;125], 'AbsTol', 1e-9)
    verifyEqual(testCase, bottlenecks.capUsage(order), [0;0;0;0;0.125], 'AbsTol', 1e-9)
    verifyEqual(testCase, bottlenecks.upperBound(order), repmat(1000, 5, 1), 'AbsTol', 1e-9)

    top3 = getEnzymeBottlenecks(ecModel, 'top', 3);
    verifyEqual(testCase, top3.uniprot, {'P1';'P2';'P3'})
end

