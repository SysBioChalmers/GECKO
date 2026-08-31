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
    % Verifies getECfromGEM copies EC codes from the conventional GEM's
    % model.eccodes into ecModel.ec.eccodes, for both the full and light
    % ecModel variants, either for every ec.rxns entry (ecRxns omitted) or
    % only for a caller-supplied logical/index subset (ecRxns given), with
    % all excluded entries left empty.
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
    % Verifies getECfromDatabase looks up EC codes from KEGG/UniProt data
    % (via loadDatabases, in 'display' mode so results are not written back
    % to model files) and stores them in ecModel.ec.eccodes, for both the
    % full and light ecModel variants, either for every ec.rxns entry
    % (ecRxns omitted) or only for a caller-supplied logical/index subset
    % (ecRxns given), with all excluded entries left empty.
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
    % Round-trips an ecModel through saveEcModel/loadEcModel and verifies the
    % reloaded model matches the original, then verifies loadConventionalGEM
    % reproduces the underlying conventional GEM, aside from the
    % annotation/date/description/version fields it does not restore and the
    % geneShortNames field it adds.
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
    % Verifies fuzzyKcatMatching against ecTestGEM's fixture BRENDA data, for
    % both the full and light ecModel variants and for both the complete
    % ec.rxns set (ecRxns omitted) and a caller-supplied subset (ecRxns
    % given). Checks the returned kcats, substrates, eccodes, wildcardLvl and
    % origin fields, confirming that a substrate match is preferred over an
    % organism match and that a wildcarded EC match is used only as a last
    % resort.
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
    % Verifies that when two ec.enzymes rows map to the same prot_<accession>
    % metabolite (distinct genes sharing one UniProt entry), applyKcatConstraints
    % sums both enzymes' protein-cost contributions into that metabolite's S
    % matrix entry, instead of only keeping the last one written.
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
    % Verifies reportEnzymeUsage's 'topAbsUsage' Name-Value argument does not
    % error when the model has fewer enzymes than requested: the default (10)
    % on ecTestGEM's 5-enzyme model, the two "report all enzymes" sentinel
    % values (0 and Inf), and an ordinary in-range value (2) must all
    % complete without indexing past the available rows.
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

    report = reportEnzymeUsage(ecModel, usageData);
    verifyEqual(testCase, height(report.topAbsUsage), 0)

    report = reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', 0);
    verifyEqual(testCase, height(report.topAbsUsage), 0)

    report = reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', Inf);
    verifyEqual(testCase, height(report.topAbsUsage), 0)

    report = reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', 2);
    verifyEqual(testCase, height(report.topAbsUsage), 0)
end


function testSigmaFitterModelMatchesReturnedSigma_tc0020(testCase)
    % Verifies that the model sigmaFitter returns has its protein-pool
    % upper bound sized to the sigma value it also returns: the
    % prot_pool_exchange upper bound must equal Ptot*f*sigma*1000 at that
    % returned sigma, not at whichever sigma the internal grid search
    % happened to try last.
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
    % Verifies that fuzzyKcatMatching's wildcard escalation matches BRENDA EC
    % codes by prefix, not by substring: escalating '1.1.1.1' to the
    % broadest wildcard level searches for EC codes starting with '1.', so an
    % unrelated code that merely contains '1.' elsewhere (e.g. '4.2.1.1')
    % must not match. The fixture pairs one genuine class-1 BRENDA entry with
    % one such unrelated, larger-kcat entry; the returned kcat must come from
    % the genuine prefix match, not the larger unrelated value.
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
    fid = fopen(fullfile(brendaDir, 'kcat.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - kcat in 1/s\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tkcat_max\tkcat_median\tn\treferences\n');
    fprintf(fid, '1.9.9.9\tm1\ttestus testus\t42\t42\t1\t*\n');
    fprintf(fid, '4.2.1.1\tm1\ttestus testus\t500\t500\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'mw.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - molecular weight in g/mol\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tmw\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'sa.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - specific activity in umol/min/mg\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tsa_max\tsa_median\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t1\t*\n');
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
    % Verifies fuzzyKcatMatching's cross-EC-token tie-break: when a
    % reaction's multiple EC tokens match BRENDA at the same origin (organism
    % tier) but different wildcard levels, the token matched with fewer
    % wildcards must win, even if a more-wildcarded token happens to carry a
    % larger kcat. The fixture gives one reaction two EC tokens that both
    % resolve at the same origin -- one an exact (0-wildcard) match with a
    % small kcat, the other a 2-wildcard match with a larger kcat -- and
    % expects the exact match's kcat and wildcardLvl to be returned.
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
    fid = fopen(fullfile(brendaDir, 'kcat.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - kcat in 1/s\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tkcat_max\tkcat_median\tn\treferences\n');
    fprintf(fid, '9.9.9.1\tunrelatedsubstrate\ttestus testus\t42\t42\t1\t*\n');
    fprintf(fid, '9.9.7.5\tunrelatedsubstrate\ttestus testus\t500\t500\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'mw.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - molecular weight in g/mol\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tmw\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'sa.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - specific activity in umol/min/mg\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tsa_max\tsa_median\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t1\t*\n');
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
    % Verifies writeOpenKineticsPredictorInput correctly maps a
    % caller-supplied ecRxns subset back to the right ec.rxns rows even when
    % that subset excludes entries that fall between the requested reactions
    % in ec.rxns order: requesting a subset whose reactions are not
    % contiguous in ec.rxns (excluding two entries in between) must still
    % write the substrate/gene-sequence pairs belonging to the requested
    % reactions, not those of the skipped-over entries.
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
    % Verifies getEnzymeBottlenecks solves the model, ranks enzymes by the
    % absolute shadow price of their prot_<id> mass-balance constraint, and
    % returns the top N rows. Checks the primal-derived columns (flux,
    % capUsage, upperBound) exactly, that rows are sorted by non-increasing
    % |shadowPrice|, and that 'top' limits the row count -- but not the exact
    % shadow-price value for an enzyme carrying zero flux at the optimum,
    % since that value is LP-dual-degenerate and can vary with the solver's
    % platform/presolve path.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    model.lb(strcmp(model.rxns,'R2')) = 0; model.ub(strcmp(model.rxns,'R2')) = 0;
    model.lb(strcmp(model.rxns,'R4')) = 0; model.ub(strcmp(model.rxns,'R4')) = 0;
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 10;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, 0.5, 0.5, 0.5, adapter);

    sol = solveLP(ecModel);
    verifyEqual(testCase, sol.stat, 1)
    verifyEqual(testCase, sol.f, 50, 'AbsTol', 1e-9)

    bottlenecks = getEnzymeBottlenecks(ecModel);
    verifyEqual(testCase, height(bottlenecks), 5)

    % Primal-derived values: unique regardless of dual degeneracy.
    [~, order] = ismember({'P1','P2','P3','P4','P5'}, bottlenecks.uniprot);
    verifyEqual(testCase, bottlenecks.flux(order), [0;0;0;55.555555555556;69.444444444444], 'AbsTol', 1e-9)
    verifyEqual(testCase, bottlenecks.capUsage(order), [0;0;0;0.055555555555556;0.069444444444444], 'AbsTol', 1e-9)
    verifyEqual(testCase, bottlenecks.upperBound(order), repmat(1000, 5, 1), 'AbsTol', 1e-9)

    % Sorted by non-increasing |shadow price| --- true whichever degenerate
    % vertex the solver lands on.
    verifyTrue(testCase, all(diff(abs(bottlenecks.shadowPrice)) <= 1e-9))

    % top limits the row count to exactly what was asked for, whatever the
    % tie-break among degenerate rows put there.
    top2 = getEnzymeBottlenecks(ecModel, 'top', 2);
    verifyEqual(testCase, height(top2), 2)
end


function testPfbaEnzymesMinimisesEnzymeUsage_tc0025(testCase)
    % Verifies pfbaEnzymes fixes the current objective (growth, via R5) as a
    % constraint and then minimises total usage_prot_* flux, reaching
    % growth=90 and enzyme_usage=125 using only usage_prot_P5 (every other
    % usage reaction at 0), within the ~1e-6 relative tolerance solveLP's own
    % minFlux=1 constraint introduces. Also checks that 'fractionOfOptimum'
    % scales both the fixed growth target and the resulting enzyme usage,
    % that 'rxnId' can override which reaction is fixed as the objective, and
    % that calling it on a gecko-light model (no usage_prot_* reactions to
    % minimise over) raises an error.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 10;
    ecModel.ec.source(:) = {'manual'};
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, 0.5, 0.5, 0.5, adapter);

    solution = pfbaEnzymes(ecModel);
    verifyEqual(testCase, solution.stat, 1)
    verifyEqual(testCase, solution.objectiveValue, 90, 'RelTol', 1e-5)
    verifyEqual(testCase, solution.enzymeUsage, 125, 'RelTol', 1e-5)

    usageIdx = find(startsWith(ecModel.rxns, 'usage_prot_'));
    usageFlux = solution.x(usageIdx);
    [~, order] = ismember(strcat('usage_prot_', {'P1','P2','P3','P4','P5'}), ecModel.rxns(usageIdx));
    verifyEqual(testCase, usageFlux(order), [0;0;0;0;125], 'AbsTol', 1e-3)

    % fractionOfOptimum halves both the fixed growth target and, on this
    % fixture, the enzyme usage needed to reach it.
    solutionHalf = pfbaEnzymes(ecModel, 'fractionOfOptimum', 0.5);
    verifyEqual(testCase, solutionHalf.objectiveValue, 45, 'RelTol', 1e-5)
    verifyEqual(testCase, solutionHalf.enzymeUsage, 62.5, 'RelTol', 1e-5)

    % rxnId overrides which reaction is fixed as the objective before
    % minimising enzyme usage, without erroring.
    solutionR3 = pfbaEnzymes(ecModel, 'rxnId', 'R3');
    verifyEqual(testCase, solutionR3.stat, 1)

    % gecko-light guard: no usage_prot_<id> machinery to minimise over.
    lightModel = makeEcModel(model, true, adapter);
    lightModel = getECfromGEM(lightModel);
    try
        pfbaEnzymes(lightModel);
        raised = false;
    catch
        raised = true;
    end
    verifyTrue(testCase, raised)
end


function ecModel = fixtureForRelaxProteomicsGreedy()
% Builds the shared fixture used by the relaxProteomicsGreedy tests below:
% ecTestGEM with uniform kcat=10 and protein pool sized via
% Ptot=0.5/f=0.5/sigma=0.5, giving unconstrained growth of 90. P5 is R5's
% sole catalyst and the true bottleneck (R2/R3's own enzymes never carry
% flux on this fixture); P4 catalyses R3, which is never on the critical
% path. Concentrations are constrained to P5=10 and P4=5, which drops growth
% to 7.2 = 90*(10/125).
geckoPath = findGECKOroot;
adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
model = getGeckoTestModel();
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
ecModel.ec.kcat(:) = 10;
ecModel.ec.source(:) = {'manual'};
ecModel = applyKcatConstraints(ecModel);
ecModel = setProtPoolSize(ecModel, 0.5, 0.5, 0.5, adapter);

ecModel.ec.concs(:) = NaN;
ecModel.ec.concs(strcmp(ecModel.ec.enzymes,'P5')) = 10;
ecModel.ec.concs(strcmp(ecModel.ec.enzymes,'P4')) = 5;
ecModel = constrainEnzConcs(ecModel);
end


function testRelaxProteomicsGreedyConverges_tc0026(testCase)
    % Verifies relaxProteomicsGreedy greedily relaxes the most-shadow-priced
    % proteomics-constrained enzyme each round until minimalGrowth is
    % reached. On this fixture it converges in exactly one step, relaxing
    % only P5 (the true bottleneck; P4's shadow price is 0 since R3 is never
    % on the critical path) straight to growth=90, and records that step's
    % trace (before=7.2, after=90, shadowPrice=-0.72).
    ecModel = fixtureForRelaxProteomicsGreedy();

    solBefore = solveLP(ecModel);
    verifyEqual(testCase, solBefore.f, 7.2, 'AbsTol', 1e-9)

    result = relaxProteomicsGreedy(ecModel, 'minimalGrowth', 50);
    verifyTrue(testCase, result.converged)
    verifyEqual(testCase, result.finalGrowth, 90, 'AbsTol', 1e-9)
    verifyEqual(testCase, fieldnames(result.relaxed), {'P5'})
    verifyEqual(testCase, result.relaxed.P5, 10, 'AbsTol', 1e-9)
    verifyEqual(testCase, numel(result.trace), 1)
    verifyEqual(testCase, result.trace(1).iteration, 0)
    verifyEqual(testCase, result.trace(1).relaxedUniprot, 'P5')
    verifyEqual(testCase, result.trace(1).growthBefore, 7.2, 'AbsTol', 1e-9)
    verifyEqual(testCase, result.trace(1).growthAfter, 90, 'AbsTol', 1e-9)
    verifyEqual(testCase, result.trace(1).shadowPrice, -0.72, 'AbsTol', 1e-9)
end


function testRelaxProteomicsGreedyExhaustsCandidates_tc0027(testCase)
    % Verifies that when minimalGrowth is unreachable, relaxProteomicsGreedy
    % relaxes every eligible enzyme (P5 then P4, in shadow-price order) and
    % returns normally with converged=false, reporting finalGrowth=90, once
    % candidates run out -- distinct from the maxIterations case below,
    % which raises an error instead.
    ecModel = fixtureForRelaxProteomicsGreedy();
    result = relaxProteomicsGreedy(ecModel, 'minimalGrowth', 1000);
    verifyFalse(testCase, result.converged)
    verifyEqual(testCase, result.finalGrowth, 90, 'AbsTol', 1e-9)
    verifyEqual(testCase, sort(fieldnames(result.relaxed)), sort({'P5';'P4'}))
end


function testRelaxProteomicsGreedyMaxIterationsRaises_tc0028(testCase)
    % Verifies that maxIterations=1 stops relaxProteomicsGreedy after
    % relaxing only P5 (growth=90, still below the unreachable target) while
    % an eligible candidate (P4) remains: this raises an error rather than
    % returning as if converged or exhausted.
    ecModel = fixtureForRelaxProteomicsGreedy();
    try
        relaxProteomicsGreedy(ecModel, 'minimalGrowth', 1000, 'maxIterations', 1);
        raised = false;
    catch
        raised = true;
    end
    verifyTrue(testCase, raised)
end


function testResolveOkpApiKeyResolutionOrder_tc0029(testCase)
    % resolveOkpApiKey: argKey beats the OKP_API_KEY environment variable,
    % which beats data/okpApiKey.txt, which beats erroring when none are set.
    % Shared by submitOpenKineticsPredictor and fetchOpenKineticsPredictor;
    % tested directly since neither caller can be run end-to-end without a
    % live API key (same reasoning runDLKcat's own Docker dependency keeps
    % it out of this suite --- external-service calls aren't unit tested
    % here, only the logic that doesn't need one).
    scratchDir = fullfile(tempname);
    mkdir(fullfile(scratchDir,'data'));
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));

    oldEnv = getenv('OKP_API_KEY');
    cleanupEnv = onCleanup(@() setenv('OKP_API_KEY', oldEnv));

    % Nothing set anywhere: errors.
    setenv('OKP_API_KEY', '');
    try
        resolveOkpApiKey('', scratchDir);
        raised = false;
    catch
        raised = true;
    end
    verifyTrue(testCase, raised)

    % File only.
    fID = fopen(fullfile(scratchDir,'data','okpApiKey.txt'), 'w');
    fprintf(fID, 'ak_from_file');
    fclose(fID);
    verifyEqual(testCase, resolveOkpApiKey('', scratchDir), 'ak_from_file')

    % Env beats file.
    setenv('OKP_API_KEY', 'ak_from_env');
    verifyEqual(testCase, resolveOkpApiKey('', scratchDir), 'ak_from_env')

    % Explicit arg beats both.
    verifyEqual(testCase, resolveOkpApiKey('ak_from_arg', scratchDir), 'ak_from_arg')
end


function testFetchOpenKineticsPredictorUseStoredMatchesReader_tc0030(testCase)
    % Verifies fetchOpenKineticsPredictor's 'useStored' path (parsing an
    % already-downloaded result file, the only path testable without a live
    % API call) delegates to readOpenKineticsPredictorOutput and produces an
    % identical kcatList struct: the 'OKP-' prefix on kcatSource, and the
    % constant source 'OpenKineticsPredictor' used by every kcat-gathering
    % function.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.metSmiles = repmat({''}, numel(ecModel.mets), 1);
    ecModel.metSmiles(strcmp(ecModel.mets,'m1c')) = {'C(C1C)O'};

    scratchDir = fullfile(tempname);
    mkdir(fullfile(scratchDir,'data'));
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    fetchAdapter = adapter;
    fetchAdapter.params.path = scratchDir;

    resultFile = fullfile(scratchDir,'data','OKP_output.csv');
    fid = fopen(resultFile,'w');
    fprintf(fid,'kcat (1/s),Source kcat,Extra Info kcat,Protein Sequence,Substrate\n');
    fprintf(fid,'12.5,Prediction from CataPro,,MRAL,C(C1C)O\n');
    fclose(fid);

    [done, viaFetch] = fetchOpenKineticsPredictor(ecModel, 'useStored', true, 'modelAdapter', fetchAdapter);
    verifyTrue(testCase, done)

    viaReader = readOpenKineticsPredictorOutput(ecModel, 'outFile', resultFile, 'modelAdapter', fetchAdapter);
    verifyEqual(testCase, viaFetch, viaReader)
    verifyEqual(testCase, viaFetch.kcatSource, {'OKP-CataPro'})
    verifyEqual(testCase, viaFetch.source, 'OpenKineticsPredictor')
end


function testGetECfromGEMRejectsMalformedComponents_tc0031(testCase)
    % Verifies getECfromGEM only accepts an EC code whose four
    % dot-separated components are each purely digits or the wildcard '-':
    % a malformed component containing letters or underscores (e.g.
    % '1.1.1.n12' or '1_2_3_4') is rejected and left empty, while a
    % well-formed wildcarded code (e.g. '1.2.3.-') is kept.
    clear model
    model.rxns = {'R1';'R2';'R3';'R4'};
    model.eccodes = {'1.1.1.1'; '1.1.1.n12'; '1_2_3_4'; '1.2.3.-'};
    model.ec.rxns = model.rxns;
    model.ec.geckoLight = false;

    ecModel = getECfromGEM(model);
    expected = {'1.1.1.1'; ''; ''; '1.2.3.-'};
    verifyEqual(testCase, ecModel.ec.eccodes, expected)
end


function testFindECInDBIntersectionDedupesWildcardPair_tc0032(testCase)
    % Verifies that reconciling EC annotations across an AND-complex's genes
    % via findECInDB de-duplicates a code that matches via both itself and
    % its own wildcard parent. Gene 1's protein carries two EC tokens,
    % '1.1.1.1' and its own wildcard parent '1.1.1.-' (space-separated, the
    % format getECstring expects for one protein's multiple activities);
    % gene 2's protein carries just '1.1.1.1'. The intersection across both
    % genes must return '1.1.1.1' once, not twice.
    geneHashMap = containers.Map({'G1','G2'}, {1,2});
    geneIndex = {1; 2};
    DBecNum = {'1.1.1.1 1.1.1.-'; '1.1.1.1'};
    DBMW = [10000; 20000];

    EC = findECInDB({'G1','G2'}, DBecNum, DBMW, geneIndex, geneHashMap);
    verifyEqual(testCase, EC, '1.1.1.1')
end


function testFuzzyKcatMatchingWildcardExhaustionDoesNotCrash_tc0033(testCase)
    % Verifies that when an EC code has no BRENDA match at any wildcard
    % level, including the fully wildcarded '-.-.-.-', fuzzyKcatMatching's
    % escalation loop terminates cleanly and reports "no match" (kcat 0,
    % origin NaN) instead of erroring.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    % An EC class unrelated, at every wildcard level, to anything in the
    % fixture BRENDA file below.
    ecModel.ec.eccodes{strcmp(ecModel.ec.rxns, 'R2_EXP_1')} = '9.9.9.9';

    scratchDir = fullfile(tempname);
    brendaDir = fullfile(scratchDir, 'data');
    mkdir(brendaDir);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    % Deliberately unrelated to '9.9.9.9' at every wildcard level, so the
    % escalation loop runs out of levels without ever matching.
    fid = fopen(fullfile(brendaDir, 'kcat.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - kcat in 1/s\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tkcat_max\tkcat_median\tn\treferences\n');
    fprintf(fid, '5.5.5.5\tm1\ttestus testus\t42\t42\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'mw.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - molecular weight in g/mol\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tmw\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(brendaDir, 'sa.tsv'), 'w');
    fprintf(fid, '# BRENDA release test-fixture - specific activity in umol/min/mg\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tsa_max\tsa_median\tn\treferences\n');
    fprintf(fid, '0.0.0.0\t*\tplaceholder\t1\t1\t1\t*\n');
    fclose(fid);
    % getPhylDistStructPath also resolves under params.path/data; reuse
    % ecTestGEM's own real fixture rather than fabricating a second one.
    copyfile(fullfile(geckoPath,'test','unit_tests','ecTestGEM','data','PhylDist.mat'), ...
        fullfile(brendaDir, 'PhylDist.mat'));

    fuzzyAdapter = adapter;
    fuzzyAdapter.params.path = scratchDir;

    % Must complete without erroring, and report a clean "no match" rather
    % than crashing or leaving a garbled value.
    kcatList = fuzzyKcatMatching(ecModel, ismember(ecModel.ec.rxns,{'R2_EXP_1'}), fuzzyAdapter);
    verifyEqual(testCase, kcatList.kcats, 0)
    verifyTrue(testCase, isnan(kcatList.origin))
end


function testGetStandardKcatUsesSubsystemKcatWheneverAnyMatches_tc0034(testCase)
    % Verifies that getStandardKcat assigns a GPR-less reaction its own
    % subsystem's kcat whenever any subsystem has a computed kcat, even when
    % the model has more than one distinct subsystem. R1 has no GPR and
    % shares subSystem 'SubA' with R3 (kcat 50); R5 sits in a different
    % subSystem 'SubB' (kcat 999, deliberately large so it would skew the
    % model-wide fallback if picked by mistake). R1 must be assigned SubA's
    % kcat (50), not the model-wide standardKcat fallback.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    ecModel.ec.kcat(:) = 0;
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns,'R3')) = 50;
    ecModel.ec.kcat(strcmp(ecModel.ec.rxns,'R5')) = 999;
    ecModel.ec.source(:) = {'manual'};

    ecModel.subSystems = repmat({{''}}, size(ecModel.rxns));
    ecModel.subSystems(strcmp(ecModel.rxns,'R1')) = {{'SubA'}};
    ecModel.subSystems(strcmp(ecModel.rxns,'R3')) = {{'SubA'}};
    ecModel.subSystems(strcmp(ecModel.rxns,'R5')) = {{'SubB'}};

    ecModel = getStandardKcat(ecModel, 'modelAdapter', adapter, 'threshold', 1, 'fillZeroKcat', false);

    r1kcat = ecModel.ec.kcat(strcmp(ecModel.ec.rxns,'R1'));
    verifyEqual(testCase, r1kcat, 50)
end


function testReadDLKcatOutputSubstrateMatchIsCaseInsensitive_tc0035(testCase)
    % Verifies readDLKcatOutput matches a DLKcat output file's substrate
    % names against model.metNames case-insensitively, so an output file
    % whose substrate names differ only in case from model.metNames (e.g.
    % 'M1' vs. 'm1') is still recognized rather than reported as
    % "substrate not found".
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);

    scratchDir = fullfile(tempname);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    mkdir(scratchDir);
    outFile = fullfile(scratchDir,'DLKcat_output_test.tsv');
    fid = fopen(outFile,'w');
    fprintf(fid,'rxns\tgenes\tsubstrates\tsmiles\tsequence\tkcat\n');
    % 'M1' (uppercase) against the model's own lowercase 'm1' -- must still match.
    fprintf(fid,'R2_EXP_1\tG1\tM1\tsmiles\tseq\t12.5\n');
    fclose(fid);

    kcatList = readDLKcatOutput(ecModel, 'outFile', outFile, 'modelAdapter', adapter);
    verifyEqual(testCase, kcatList.kcats, 12.5)
    verifyEqual(testCase, kcatList.substrates, {'M1'})
end


function testMakeEcModelSortsEcGenes_tc0036(testCase)
    % Verifies makeEcModel sorts ec.genes alphabetically (permuting
    % ec.enzymes/ec.mw/ec.sequence in lockstep), regardless of the order
    % genes are declared in model.genes. Reorders ecTestGEM's genes (which
    % are otherwise already alphabetical) before building the ecModel, so
    % the test actually exercises the sort rather than an already-sorted
    % pass-through.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    newOrder = [5 3 1 4 2];
    model.genes = model.genes(newOrder);
    model.rxnGeneMat = model.rxnGeneMat(:,newOrder);

    ecModel = makeEcModel(model, false, adapter);

    expEcGenes = {'G1';'G2';'G3';'G4';'G5'};
    verifyEqual(testCase,ecModel.ec.genes,expEcGenes)
    expEnzymes = {'P1';'P2';'P3';'P4';'P5'};
    verifyEqual(testCase,ecModel.ec.enzymes,expEnzymes)
    expMW = [10000;20000;30000;40000;50000];
    verifyEqual(testCase,ecModel.ec.mw,expMW)
    expSequence = {'MRAL';'MNTD';'MSYN';'MDFM';'MLFK'};
    verifyEqual(testCase,ecModel.ec.sequence,expSequence);
end


function testReadDLKcatOutputUnrecognizedSubstrateErrorsCleanly_tc0037(testCase)
    % Verifies that readDLKcatOutput reports a genuinely unrecognized
    % substrate with a clean error naming the substrate, rather than
    % crashing.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);

    scratchDir = fullfile(tempname);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    mkdir(scratchDir);
    outFile = fullfile(scratchDir,'DLKcat_output_test.tsv');
    fid = fopen(outFile,'w');
    fprintf(fid,'rxns\tgenes\tsubstrates\tsmiles\tsequence\tkcat\n');
    fprintf(fid,'R2_EXP_1\tG1\tNOTAREALMETABOLITE\tsmiles\tseq\t12.5\n');
    fclose(fid);

    try
        readDLKcatOutput(ecModel, 'outFile', outFile, 'modelAdapter', adapter);
        errored = false;
        msg = '';
    catch ME
        errored = true;
        msg = ME.message;
    end
    verifyTrue(testCase, errored)
    verifyFalse(testCase, contains(msg, 'dispEM'))
    verifyTrue(testCase, contains(msg, 'NOTAREALMETABOLITE'))
end


function testLoadDatabasesDuplicateUniprotEntriesErrorsCleanly_tc0038(testCase)
    % Verifies that loadDatabases reports duplicate UniProt entries with a
    % clean error, rather than crashing.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));

    scratchDir = fullfile(tempname);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    mkdir(fullfile(scratchDir,'data'));
    fid = fopen(fullfile(scratchDir,'data','uniprot.tsv'),'w');
    fprintf(fid,'Entry\tGene Names (ordered locus)\tEC number\tMass\tSequence\n');
    fprintf(fid,'P1\tG1\t1.1.1.1\t10000\tMRAL\n');
    fprintf(fid,'P1\tG1\t1.1.1.1\t10000\tMRAL\n');
    fclose(fid);

    dbAdapter = adapter;
    dbAdapter.params.path = scratchDir;

    try
        loadDatabases('uniprot', dbAdapter);
        errored = false;
        msg = '';
    catch ME
        errored = true;
        msg = ME.message;
    end
    verifyTrue(testCase, errored)
    verifyFalse(testCase, contains(msg, 'dispEM'))
    verifyTrue(testCase, contains(msg, 'Duplicate entries'))
end


function testAddNewRxnsToECMissingGeneErrorsCleanly_tc0039(testCase)
    % Verifies that addNewRxnsToEC reports a gene referenced in grRules but
    % missing from newEnzymes with a clean error naming that gene, rather
    % than crashing.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    ecModel = makeEcModel(model, false, adapter);

    newRxns.rxns      = {'RNEWA'};
    newRxns.rxnNames  = {'RNEWA'};
    newRxns.equations = {'m1[c] => e2[e]'};
    newRxns.grRules   = {'GNEW1 or GNEW2'};

    % GNEW2 is referenced in grRules but deliberately left out of newEnzymes below.
    newEnzymes.enzymes = {'ENEW1'};
    newEnzymes.genes   = {'GNEW1'};
    newEnzymes.mw      = 15000;

    try
        addNewRxnsToEC(ecModel, newRxns, newEnzymes, adapter);
        errored = false;
        msg = '';
    catch ME
        errored = true;
        msg = ME.message;
    end
    verifyTrue(testCase, errored)
    verifyFalse(testCase, contains(msg, 'dispEM'))
    verifyTrue(testCase, contains(msg, 'GNEW2'))
end


function testGetSubsetEcModelMismatchedReactionsErrorsCleanly_tc0040(testCase)
    % Verifies that getSubsetEcModel reports reactions in smallGEM that
    % cannot be matched to bigEcModel with a clean error naming those
    % reactions, rather than crashing.
    clear bigEcModel smallGEM
    bigEcModel.rxns = {'R1';'R2_EXP_1'};
    smallGEM.rxns   = {'R1';'R3'};

    try
        getSubsetEcModel(bigEcModel, smallGEM);
        errored = false;
        msg = '';
    catch ME
        errored = true;
        msg = ME.message;
    end
    verifyTrue(testCase, errored)
    verifyFalse(testCase, contains(msg, 'dispEM'))
    verifyTrue(testCase, contains(msg, 'R3'))
end


function testMergeKcatsGeneralizesToNSourcesWithExactTier_tc0041(testCase)
    % mergeKcats generalizes mergeDLKcatAndFuzzyKcats (kept as a thin,
    % deprecated wrapper around it) from two fixed inputs to any number of
    % kcatLists, each of which may itself mix several sources via a per-row
    % kcatSource field -- e.g. fetchOpenKineticsPredictor's own output,
    % modelled here as okpList, carrying both an exact BRENDA hit (no
    % wildcard/origin metadata at all) and a CataPro prediction.
    fuzzyList.rxns        = {'R1'; 'R2'};
    fuzzyList.kcats       = [5; 8];
    fuzzyList.genes       = {{}; {}};
    fuzzyList.substrates  = {{}; {}};
    fuzzyList.eccodes     = {'1.1.1.1'; '1.1.1.2'};
    fuzzyList.wildcardLvl = [0; 0];
    fuzzyList.origin      = [1; 1];
    fuzzyList.source      = 'brenda';

    dlkcatList.rxns       = {'R1'; 'R3'};
    dlkcatList.kcats      = [50; 30];
    dlkcatList.genes      = {'G1'; 'G3'};
    dlkcatList.substrates = {'m1'; 'm3'};
    dlkcatList.source     = 'DLKcat';

    okpList.rxns       = {'R1'; 'R4'};
    okpList.kcats      = [99; 15];
    okpList.genes      = {[]; []};
    okpList.substrates = {[]; []};
    okpList.kcatSource = {'BRENDA'; 'CataPro'};

    mergedKcatList = mergeKcats({fuzzyList, dlkcatList, okpList}, ...
        {'database_exact', 'database_top', 'dlkcat'});

    % R1: fuzzy gives 'database_top' (wc=0, origin=1); OKP's BRENDA row has
    % no wildcard/origin metadata at all, so it is an exact hit --
    % 'database_exact' outranks 'database_top', so OKP's 99 wins, not
    % fuzzy's 5 or DLKcat's 50 (DLKcat is never even reached for R1).
    r1 = strcmp(mergedKcatList.rxns, 'R1');
    verifyEqual(testCase, sum(r1), 1)
    verifyEqual(testCase, mergedKcatList.kcats(r1), 99)

    % R2: only in fuzzy ('database_top') -- kept.
    r2 = strcmp(mergedKcatList.rxns, 'R2');
    verifyEqual(testCase, mergedKcatList.kcats(r2), 8)

    % R3: only in DLKcat -- kept via the 'dlkcat' tier.
    r3 = strcmp(mergedKcatList.rxns, 'R3');
    verifyEqual(testCase, mergedKcatList.kcats(r3), 30)

    % R4: only in OKP as 'CataPro', a source not listed in sourcePriority --
    % dropped, not present in the merged output at all.
    verifyEqual(testCase, sum(strcmp(mergedKcatList.rxns, 'R4')), 0)
end


function testMergeDLKcatAndFuzzyKcatsDelegatesToMergeKcats_tc0042(testCase)
    % Verifies mergeDLKcatAndFuzzyKcats, called positionally, delegates to
    % mergeKcats using a fixed three-tier priority order (database_top >
    % dlkcat > database_bottom), using the same style of fixture
    % testKcats_tc0011 already established for this function.
    fuzzyList.rxns        = {'R1'; 'R2'; 'R3'};
    fuzzyList.kcats       = [5; 8; 2];
    fuzzyList.genes       = {{}; {}; {}};
    fuzzyList.substrates  = {{}; {}; {}};
    fuzzyList.eccodes     = {'1.1.1.1'; '1.1.1.2'; '1.1.1.3'};
    % R1: origin 1, top-tier match. R2: no BRENDA row at all (DLKcat only).
    % R3: wildcardLvl 1, only within the bottom tier's wildcard allowance.
    fuzzyList.wildcardLvl = [0; NaN; 1];
    fuzzyList.origin      = [1; NaN; 6];
    fuzzyList.source      = 'brenda';

    dlkcatList.rxns       = {'R1'; 'R2'};
    dlkcatList.kcats      = [50; 80];
    dlkcatList.genes      = {'G1'; 'G2'};
    dlkcatList.substrates = {'m1'; 'm2'};
    dlkcatList.source     = 'DLKcat';

    mergedKcatList = mergeDLKcatAndFuzzyKcats(dlkcatList, fuzzyList);

    % R1: BRENDA top-tier match (5) wins over DLKcat (50).
    r1 = strcmp(mergedKcatList.rxns, 'R1');
    verifyEqual(testCase, mergedKcatList.kcats(r1), 5)
    % R2: no BRENDA row -- DLKcat (80) is the only candidate.
    r2 = strcmp(mergedKcatList.rxns, 'R2');
    verifyEqual(testCase, mergedKcatList.kcats(r2), 80)
    % R3: BRENDA bottom-tier wildcard match (2), no DLKcat row to compete.
    r3 = strcmp(mergedKcatList.rxns, 'R3');
    verifyEqual(testCase, mergedKcatList.kcats(r3), 2)
end


function testLoadBRENDAdataParsesNewTsvFormat_tc0043(testCase)
    % Verifies loadBRENDAdata parses kcat.tsv/mw.tsv/sa.tsv files (a
    % `#`-prefixed release line, a tab-delimited column header, bare EC
    % codes, plain organism strings, and, for kcat/sa, both a _max and
    % _median aggregate per row), and that the _max column -- not _median --
    % is the value used.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));

    scratchDir = fullfile(tempname);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    mkdir(fullfile(scratchDir,'data'));
    fid = fopen(fullfile(scratchDir,'data','kcat.tsv'), 'w');
    fprintf(fid, '# BRENDA release 2026.1 generated 2026-05-28 - CC BY 4.0 - kcat in 1/s\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tkcat_max\tkcat_median\tn\treferences\n');
    fprintf(fid, '1.2.3.4\tm1\ttestorg\t99\t5\t3\tPMID:1\n');
    fclose(fid);
    fid = fopen(fullfile(scratchDir,'data','mw.tsv'), 'w');
    fprintf(fid, '# BRENDA release 2026.1 generated 2026-05-28 - CC BY 4.0 - molecular weight in g/mol\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tmw\tn\treferences\n');
    fprintf(fid, '1.2.3.4\t*\ttestorg\t60000\t1\t*\n');
    fclose(fid);
    fid = fopen(fullfile(scratchDir,'data','sa.tsv'), 'w');
    fprintf(fid, '# BRENDA release 2026.1 generated 2026-05-28 - CC BY 4.0 - specific activity in umol/min/mg\n');
    fprintf(fid, 'ec_code\tsubstrate\torganism\tsa_max\tsa_median\tn\treferences\n');
    fprintf(fid, '1.2.3.4\t*\ttestorg\t120\t7\t2\tPMID:2\n');
    fclose(fid);

    brendaAdapter = adapter;
    brendaAdapter.params.path = scratchDir;

    [KCATcell, SAcell] = loadBRENDAdata('modelAdapter', brendaAdapter);

    verifyEqual(testCase, KCATcell{1}, {'1.2.3.4'})
    verifyEqual(testCase, KCATcell{2}, {'m1'})
    verifyEqual(testCase, KCATcell{3}, {'testorg'})
    verifyEqual(testCase, KCATcell{4}, 99, 'AbsTol', 1e-9) % kcat_max, not kcat_median (5)

    verifyEqual(testCase, SAcell{1}, {'1.2.3.4'})
    verifyEqual(testCase, SAcell{2}, {'testorg'})
    % sa_max [umol/min/mg -> mmol/s/g] (x1/60) times mw [g/mol -> g/mmol] (x1/1000)
    expectedKcat = (120/60) * (60000/1000);
    verifyEqual(testCase, SAcell{3}, expectedKcat, 'AbsTol', 1e-9)
    verifyEqual(testCase, SAcell{4}, 60000/1000, 'AbsTol', 1e-9)
end


function testMakeEcModelFallsBackToKeggForUnmatchedGenes_tc0044(testCase)
    % Verifies that makeEcModel falls back to databases.kegg for a gene with
    % no UniProt match: ec.mw and ec.sequence are filled from the KEGG row,
    % and ec.enzymes gets the UniProt accession carried on that row, or --
    % when that column is itself empty, as tested here -- the bare KEGG gene
    % id as a fallback identifier. Such a gene must not appear in the
    % noUniprot output.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();

    scratchDir = fullfile(tempname);
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    mkdir(fullfile(scratchDir,'data'));
    % G5 is deliberately absent from uniprot.tsv (only G1-G4 are UniProt-matched)
    % but present in kegg.tsv with an empty UniProt-accession column (col 1),
    % forcing the bare-KEGG-gene-id fallback rather than the more common
    % KEGG-row-carries-a-UniProt-accession case.
    fid = fopen(fullfile(scratchDir,'data','uniprot.tsv'), 'w');
    fprintf(fid, 'Entry\tGene Names (ordered locus)\tEC number\tMass\tSequence\n');
    fprintf(fid, 'P1\tG1\t1.1.1.1\t10000\tMRAL\n');
    fprintf(fid, 'P2\tG2\t1.1.1.1\t20000\tMNTD\n');
    fprintf(fid, 'P3\tG3\t1.1.1.1\t30000\tMSYN\n');
    fprintf(fid, 'P4\tG4\t1.1.2.1\t40000\tMDFM\n');
    fclose(fid);
    fid = fopen(fullfile(scratchDir,'data','kegg.tsv'), 'w');
    fprintf(fid, ',G5,G5,1.1.1.3,50000,"",MLFK\n');
    fclose(fid);

    keggAdapter = adapter;
    keggAdapter.params.path = scratchDir;

    [ecModel, noUniprot] = makeEcModel(model, false, keggAdapter);

    g5 = strcmp(ecModel.ec.genes,'G5');
    verifyTrue(testCase, any(g5))
    verifyEqual(testCase, ecModel.ec.enzymes{g5}, 'G5') % bare KEGG gene id, no UniProt accession on the row
    verifyEqual(testCase, ecModel.ec.mw(g5), 50000)
    verifyEqual(testCase, ecModel.ec.sequence{g5}, 'MLFK')
    % G5 was resolved (via KEGG), so it must not appear in noUniprot.
    verifyEqual(testCase, sum(strcmp(noUniprot,'G5')), 0)
end


function testReportEnzymeUsageSkipsOnlyFullyInactiveEnzymes_tc0045(testCase)
    % Verifies reportEnzymeUsage's topAbsUsage table skips an enzyme only
    % when every one of its reactions carries no flux; an enzyme with some
    % flux-carrying and some inactive reactions is still reported, attributed
    % only to the flux-carrying subset. Uses a minimal hand-built
    % ecModel-like struct with three enzymes covering all three cases:
    %   P1 -- Ra, Rb, neither carries flux            -> skipped entirely
    %   P2 -- Rc, Rd; only Rc carries flux             -> single-branch, Rc only
    %   P3 -- Re, Rf, Rg; Re and Rf carry flux, Rg not -> combined-branch, Re+Rf only
    ecModel.rxns     = {'Ra';'Rb';'Rc';'Rd';'Re';'Rf';'Rg';'prot_pool_exchange'};
    ecModel.rxnNames = ecModel.rxns;
    ecModel.grRules  = {'';'';'';'';'';'';'';''};
    ecModel.mets     = {'prot_P3'};
    ecModel.S        = zeros(1,numel(ecModel.rxns));
    ecModel.S(1,strcmp(ecModel.rxns,'Re')) = -1;
    ecModel.S(1,strcmp(ecModel.rxns,'Rf')) = -1;
    ecModel.ub       = 1000*ones(numel(ecModel.rxns),1);
    ecModel.ub(strcmp(ecModel.rxns,'prot_pool_exchange')) = 100;

    ecModel.ec.enzymes = {'P1';'P2';'P3'};
    ecModel.ec.genes   = {'G1';'G2';'G3'};
    ecModel.ec.rxns    = {'Ra';'Rb';'Rc';'Rd';'Re';'Rf';'Rg'};
    ecModel.ec.kcat    = [1;1;2;2;3;3;3];
    ecModel.ec.source  = repmat({'manual'},7,1);
    ecModel.ec.rxnEnzMat = zeros(7,3);
    ecModel.ec.rxnEnzMat(1,1) = 1; % Ra -> P1
    ecModel.ec.rxnEnzMat(2,1) = 1; % Rb -> P1
    ecModel.ec.rxnEnzMat(3,2) = 1; % Rc -> P2
    ecModel.ec.rxnEnzMat(4,2) = 1; % Rd -> P2
    ecModel.ec.rxnEnzMat(5,3) = 1; % Re -> P3
    ecModel.ec.rxnEnzMat(6,3) = 1; % Rf -> P3
    ecModel.ec.rxnEnzMat(7,3) = 1; % Rg -> P3

    usageData.protID   = {'P1';'P2';'P3'};
    usageData.absUsage = [0;4;5];
    usageData.capUsage = [0;0.4;0.5];
    fluxes = zeros(numel(ecModel.rxns),1);
    fluxes(strcmp(ecModel.rxns,'Rc')) = 4;
    fluxes(strcmp(ecModel.rxns,'Re')) = 2;
    fluxes(strcmp(ecModel.rxns,'Rf')) = 3;
    usageData.fluxes = fluxes;

    report = reportEnzymeUsage(ecModel, usageData);
    top = report.topAbsUsage;

    % P1: fully inactive -- must not appear at all, not even as a placeholder.
    verifyEqual(testCase, sum(strcmp(top.protID,'P1')), 0)
    % P2 contributes 1 row (single-branch); P3 contributes 3 (a '===' header
    % row plus one detail row per active reaction, Re and Rf -- the combined
    % branch's pre-existing, unrelated behaviour for any enzyme with more
    % than one flux-carrying reaction).
    verifyEqual(testCase, height(top), 4)

    % P2: mixed, one active reaction -- reported via the single-branch,
    % attributed to Rc (not a placeholder, not Rd).
    p2 = strcmp(top.protID,'P2');
    verifyEqual(testCase, sum(p2), 1)
    verifyEqual(testCase, top.rxnID{p2}, 'Rc')
    verifyEqual(testCase, top.kcat(p2), 2)
    verifyEqual(testCase, top.absUsage(p2), 4)

    % P3: mixed, two active reactions out of three -- reported via the
    % combined-branch (header row + one detail row per active reaction),
    % split only across Re and Rf (not Rg).
    p3 = strcmp(top.protID,'P3');
    verifyEqual(testCase, sum(p3), 3)
    p3header = p3 & strcmp(top.rxnID,'===');
    p3detail = p3 & ~strcmp(top.rxnID,'===');
    verifyEqual(testCase, sum(p3header), 1)
    verifyEqual(testCase, top.absUsage(p3header), 5, 'AbsTol', 1e-9) % enzyme-level total
    verifyEqual(testCase, sort(top.rxnID(p3detail)), sort({'Re';'Rf'}))
    verifyEqual(testCase, sum(top.absUsage(p3detail)), 5, 'AbsTol', 1e-9) % detail rows split it
end


function testCalculateMWReconciledWithGeckopy_tc0046(testCase)
    % Verifies calculateMW's residue-mass handling: the water constant is
    % 18.01528; ambiguous residues X and B are computed as the mean of their
    % possible standard residues rather than a fixed constant; residue
    % letters are matched case-insensitively; and a sequence with no
    % recognized residue returns NaN rather than just the water mass.
    verifyTrue(testCase, isnan(calculateMW('')))
    verifyTrue(testCase, isnan(calculateMW('123 ---')))
    % X: mean of the 20 standard residues (118.885) + water.
    verifyEqual(testCase, calculateMW('X'), 136.90028, 'AbsTol', 1e-9)
    % B: mean(D,N) = 114.595 + water.
    verifyEqual(testCase, calculateMW('B'), 132.61028, 'AbsTol', 1e-9)
    % Case-insensitive: all four letters count as 'A', not just the two
    % uppercase ones.
    verifyEqual(testCase, calculateMW('AaAa'), 302.33528, 'AbsTol', 1e-9)
    % The 20 standard residues once each, plus the water constant only.
    verifyEqual(testCase, calculateMW('ACDEFGHIKLMNPQRSTVWY'), 2395.71528, 'AbsTol', 1e-9)
end

function testApplyKcatConstraintsLightSkipsUnassignedIsozyme_tc0047(testCase)
    % raven-gecko-parity#57: in the light formulation, when one isozyme of
    % a reaction has no kcat assigned (kcat==0, i.e. an Inf MW/kcat) and
    % another isozyme does, applyKcatConstraints used to correct the
    % unassigned isozyme's Inf cost to 0 *before* taking min() across
    % isozymes, so that fabricated zero always won -- silently leaving
    % the reaction unconstrained by enzyme usage even though a real kcat
    % was available on another isozyme. Fixed to drop isozymes with no
    % usable kcat before selecting the cheapest one, matching geckopy's
    % apply_kcat_constraints.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    lecModel = makeEcModel(model, true, adapter);
    lecModel = getECfromGEM(lecModel);
    lecModel = applyComplexData(lecModel, [], adapter, false);

    % R2's two isozymes: 001_R2 (complex G1+G2), 002_R2 (single G3, MW 30000).
    complexIdx = find(strcmp(lecModel.ec.rxns,'001_R2'));
    singleIdx  = find(strcmp(lecModel.ec.rxns,'002_R2'));
    lecModel.ec.kcat(complexIdx) = 0; % unassigned
    lecModel.ec.kcat(singleIdx)  = 7; % real kcat

    lecModel = applyKcatConstraints(lecModel);

    % Only the single-gene isozyme has a usable kcat, so its cost must be
    % the one written -- not the fabricated zero from the unassigned
    % complex.
    expected = -30000/7/3600;
    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2'))),expected,"AbsTol",10^-10)
end

function testApplyKcatConstraintsLightNoValidIsozymeUncosted_tc0048(testCase)
    % Companion to tc0047: when *every* isozyme of a reaction lacks a
    % usable kcat, the reaction is left unconstrained (S coefficient 0),
    % same as before the raven-gecko-parity#57 fix.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    lecModel = makeEcModel(model, true, adapter);
    lecModel = getECfromGEM(lecModel);
    lecModel = applyComplexData(lecModel, [], adapter, false);

    lecModel = applyKcatConstraints(lecModel); % default kcat==0 everywhere

    verifyEqual(testCase,full(lecModel.S(strcmp(lecModel.mets, 'prot_pool'),strcmp(lecModel.rxns, 'R2'))),0,"AbsTol",10^-10)
end

function testCopyECtoGEMOverwriteNeverErasesWithEmptySource_tc0049(testCase)
    % raven-gecko-parity#79: overwrite=true used to copy ec.eccodes
    % verbatim, including empty entries -- so a reaction with a real,
    % populated eccodes annotation could be silently erased just because
    % the ec-structure's entry for the same reaction happened to be
    % empty. Fixed so an empty ec.eccodes source is never copied,
    % matching geckopy's copy_ec_to_gem.
    ecModel = struct();
    ecModel.rxns = {'R1';'R2';'R3'};
    ecModel.eccodes = {'1.1.1.1';'';'2.2.2.2'}; % R1, R3 already populated
    ecModel.ec = struct();
    ecModel.ec.rxns = {'R1';'R2';'R3'};
    ecModel.ec.eccodes = {'';'3.3.3.3';''}; % R1, R3 have no new info; R2 does

    result = copyECtoGEM(ecModel, 'overwrite', true);
    % R1's and R3's real annotations must survive an empty ec.eccodes source.
    verifyEqual(testCase,result.eccodes{1},'1.1.1.1')
    verifyEqual(testCase,result.eccodes{3},'2.2.2.2')
    % R2 had nothing before and gets the new, non-empty value.
    verifyEqual(testCase,result.eccodes{2},'3.3.3.3')
end

function testCopyECtoGEMOverwriteFalseStillFillsEmptyEntries_tc0050(testCase)
    % Sanity check that the overwrite=false path (unchanged by
    % raven-gecko-parity#79) still fills only empty entries.
    ecModel = struct();
    ecModel.rxns = {'R1';'R2'};
    ecModel.eccodes = {'1.1.1.1';''};
    ecModel.ec = struct();
    ecModel.ec.rxns = {'R1';'R2'};
    ecModel.ec.eccodes = {'9.9.9.9';'3.3.3.3'};

    result = copyECtoGEM(ecModel); % overwrite defaults to false
    verifyEqual(testCase,result.eccodes{1},'1.1.1.1') % untouched
    verifyEqual(testCase,result.eccodes{2},'3.3.3.3') % filled in
end

function testMakeEcModelSubSystemsMatchModelConvention_tc0051(testCase)
    % GECKO#412: usage_prot_* and prot_pool_exchange always assigned a
    % bare-char subSystems entry ('Protein usage'), regardless of whether
    % the rest of the model used that convention or the nested
    % cell-array-of-cell-arrays one (the common case for a real GEM).
    % Reported against a RAVEN release whose addRxns did not reconcile a
    % resulting bare-char/nested mix, so checkModelStruct rejected it at
    % save time ("the subSystems field must be a cell array"). Current
    % RAVEN's addRxns reconciles that mismatch on its own, so this test
    % cannot reproduce the save-time crash directly; it instead pins down
    % that makeEcModel itself now always emits a value matching the
    % model's own convention, rather than relying on addRxns to paper
    % over a mismatch.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();
    model.subSystems = repmat({{''}}, size(model.rxns));
    model.subSystems(strcmp(model.rxns,'R2')) = {{'Metabolism'}};

    ecModel = makeEcModel(model, false, adapter);

    isNested  = cellfun(@iscell, ecModel.subSystems);
    isCellStr = cellfun(@ischar, ecModel.subSystems);
    verifyTrue(testCase, all(isNested) || all(isCellStr), ...
        'ecModel.subSystems mixes bare chars with nested cells')

    usageIdx = find(startsWith(ecModel.rxns,'usage_prot_'));
    verifyGreaterThan(testCase, numel(usageIdx), 0)
    for i = 1:numel(usageIdx)
        verifyEqual(testCase, ecModel.subSystems{usageIdx(i)}, {'Protein usage'})
    end

    poolIdx = find(strcmp(ecModel.rxns,'prot_pool_exchange'));
    verifyEqual(testCase, ecModel.subSystems{poolIdx}, {'Protein usage'})
end

function testEcFVAIsozymeSplitReactionUsesExactDiagonalBound_tc0052(testCase)
    % raven-gecko-parity#72: for a canonical reaction split across isozyme
    % variants that compete for a shared limiting resource, ecFVA used to
    % read each variant's row-wise max/min across *every* canonical
    % group's solve before summing forward and subtracting reverse
    % variants -- an "envelope" that can combine values from different,
    % not-necessarily-jointly-feasible flux distributions. Confirmed
    % empirically here, not just argued from code: R2's two isozymes
    % (G1+G2 complex, G3 single) share one upstream supply (R1, capped at
    % 12, the only route to R2 since R3/R4 are blocked); each isozyme can
    % reach 12 alone, but never jointly past 12. The unfixed ecFVA reports
    % maxFlux(R2) = 24 -- double what any feasible flux distribution can
    % attain -- because some *other* canonical group's solve happens to
    % report R2_EXP_1 near its own max as a non-objective side effect, and
    % a different group's solve does the same for R2_EXP_2, and the
    % row-wise reduction sums the two independently. The fixed version
    % reads each canonical reaction's bound only from its own group's
    % solve (the "diagonal"), matching geckopy's ec_fva.
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    model = getGeckoTestModel();

    % Force all m1->m2 flux through R2's two isozymes.
    model.lb(strcmp(model.rxns,'R3')) = 0; model.ub(strcmp(model.rxns,'R3')) = 0;
    model.lb(strcmp(model.rxns,'R4')) = 0; model.ub(strcmp(model.rxns,'R4')) = 0;
    % Shared upstream supply, tightly capped: the two isozymes must
    % compete for it.
    model.lb(strcmp(model.rxns,'R1')) = 0; model.ub(strcmp(model.rxns,'R1')) = 12;
    % Forward-only: R2's reverse direction would let a futile
    % R2_EXP_1/R2_REV_EXP_1 cycle inflate R2_EXP_1's own flux with no net
    % mass-balance effect, unrelated to what this test is exercising.
    model.lb(strcmp(model.rxns,'R2')) = 0;

    ecModel = makeEcModel(model, false, adapter);
    ecModel = getECfromGEM(ecModel);
    % Generous kcat everywhere: the shared R1 supply must be the binding
    % constraint, not enzyme cost.
    ecModel.ec.kcat(:) = 1000;
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel, [], [], [], adapter);
    ecModel.ub(strcmp(ecModel.rxns,'prot_pool_exchange')) = 1e6;

    % ecFVA unconditionally calls parpool/gcp, which errors outright
    % (rather than just running serially) where Parallel Computing
    % Toolbox isn't available -- a pre-existing property of ecFVA.m, not
    % something this fix introduces or could work around here. Caught
    % here, rather than pre-checked with license(), because a pre-check
    % proved unreliable in practice (observed passing standalone and in
    % this same suite's own environment probe, yet still reporting
    % unavailable when this test itself ran) -- reacting to the actual
    % failure is the one check that cannot give a false negative.
    try
        [~, maxFlux] = ecFVA(ecModel, model);
    catch ME
        if strcmp(ME.identifier,'MATLAB:UndefinedFunction') && ...
                (contains(ME.message,'gcp') || contains(ME.message,'parpool'))
            testCase.assumeFail('Parallel Computing Toolbox not available; ecFVA requires it.')
        end
        rethrow(ME)
    end
    r2 = strcmp(model.rxns, 'R2');
    verifyEqual(testCase, maxFlux(r2), 12, 'AbsTol', 1e-6)
end

