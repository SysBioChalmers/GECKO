%run this test case with the command
%results = runtests('kcatAssignmentTests.m')

function tests = kcatAssignmentTests
    tests = functiontests(localfunctions);
end

function testAssignKcatsReproducesChain_tc0079(testCase)
    % The gecko_classic preset gives the same model as the steps it
    % replaces: fuzzyKcatMatching, readDLKcatOutput, mergeKcats,
    % assignKcatValues, applyCustomKcats, getKcatAcrossIsozymes,
    % getStandardKcat and applyKcatConstraints.
    [base, adapter, project] = buildProject();
    fuzzy = fuzzyKcatMatching(base, [], adapter);
    dlkcat = readDLKcatOutput(base, 'modelAdapter', adapter);
    merged = mergeKcats({fuzzy, dlkcat}, {'database_top','dlkcat','database_bottom'});
    A = assignKcatValues(base, merged);
    A = applyCustomKcats(A, 'modelAdapter', adapter);
    A = getKcatAcrossIsozymes(A);
    A = getStandardKcat(A, 'modelAdapter', adapter);
    A = applyKcatConstraints(A);

    kcatsFromBrenda(base, 'modelAdapter', adapter);
    kcatsFromDLKcat(base, 'modelAdapter', adapter);
    [B, decisions, report] = assignKcats(base, 'modelAdapter', adapter);
    verifyEqual(testCase, B.ec.rxns, A.ec.rxns)
    verifyEqual(testCase, B.ec.kcat, A.ec.kcat)
    verifyEqual(testCase, B.ec.source, A.ec.source)
    verifyEqual(testCase, B.S, A.S)
    verifyEqual(testCase, report.preset, 'gecko_classic')
    verifyEqual(testCase, decisions.ecRxn, B.ec.rxns)
    verifyTrue(testCase, isfile(fullfile(project, 'output', [char(base.id) '_kcatDecisions.tsv'])))
end

function testCollectKcatsLeavesModel_tc0080(testCase)
    % collectKcats decides without changing the model; a missing source is
    % reported; custom.tsv decides first.
    [base, adapter, project] = buildProject();
    kcatsFromBrenda(base, 'modelAdapter', adapter);
    [decisions, report, decided] = collectKcats(base, 'modelAdapter', adapter, 'writeReport', false);
    verifyEqual(testCase, base.ec.kcat, zeros(numel(base.ec.rxns), 1))
    verifyEqual(testCase, report.missing, {'dlkcat'})
    verifyTrue(testCase, any(decided.ec.kcat > 0))

    [baseRxns, enz] = getEcRxnKeys(base);
    i = find(strcmp(base.ec.rxns, 'R3'));
    writeLines(fullfile(project, 'data', 'kcats', 'custom.tsv'), {'# source: custom', '# kind: curated', ...
        'rxn	enzymes	kcat	source', sprintf('%s	%s	777	custom', baseRxns{i}, enz{i})});
    decisions = collectKcats(base, 'modelAdapter', adapter, 'writeReport', false);
    verifyEqual(testCase, decisions.kcat(i), 777)
    verifyEqual(testCase, decisions.source{i}, 'custom')
    verifyEqual(testCase, decisions.stage{i}, '1')
end

function testStagesFromToml_tc0081(testCase)
    % Stages from model_adapter.toml, and their errors.
    [base, adapter, project] = buildProject();
    kcatsFromBrenda(base, 'modelAdapter', adapter);
    kcatsFromDLKcat(base, 'modelAdapter', adapter);
    adapter.params.kcat_assignment = struct('stage', struct('action', 'select', ...
        'sources', {{'dlkcat'}}, 'criteria', 'min'), 'isozymes', false, 'standard', false);
    decisions = collectKcats(base, 'modelAdapter', adapter, 'writeReport', false);
    verifyEqual(testCase, unique(decisions.source(decisions.kcat > 0)), {'DLKcat'})
    verifyEqual(testCase, decisions.kcat(strcmp(decisions.ecRxn, 'R2_EXP_1')), 6.5)
    adapter.params.kcat_assignment.preset = 'gecko_median';
    verifyError(testCase, @() collectKcats(base, 'modelAdapter', adapter), 'kcatAssignment:config')
    adapter.params.kcat_assignment = struct('preset', 'gecko_fast');
    verifyError(testCase, @() collectKcats(base, 'modelAdapter', adapter), 'kcatAssignment:config')

    % The same stages read from a TOML file.
    toml = fullfile(project, 'model_adapter.toml');
    writeLines(toml, {'conv_gem = "models/testModel.xml"', 'org_name = "testus testus"', ...
        '[kcat_assignment]', 'isozymes = false', '[[kcat_assignment.stage]]', 'sources = ["dlkcat"]', ...
        'criteria = "min"', '[kcat_sources.brenda]', 'aggregation = "median"'});
    params = loadAdapterParams(toml);
    verifyEqual(testCase, params.kcat_assignment.stage.sources, {'dlkcat'})
    verifyEqual(testCase, params.kcat_assignment.standard, true)
    verifyEqual(testCase, params.kcat_sources.brenda.aggregation, 'median')
end

function testUpdateKcatSources_tc0082(testCase)
    % Missing and outdated sources are gathered; current ones are not.
    [base, adapter, project] = buildProject();
    adapter.params.kcat_assignment = struct('preset', 'gecko_median');
    report = updateKcatSources(base, 'dryRun', true, 'modelAdapter', adapter);
    verifyEqual(testCase, report.source', {'custom','brenda','dlkcat'})
    verifyEqual(testCase, report.status', {'skipped','missing','missing'})
    verifyFalse(testCase, isfolder(fullfile(project, 'data', 'kcats')))
    report = updateKcatSources(base, 'modelAdapter', adapter);
    verifyEqual(testCase, report.status', {'skipped','gathered','gathered'})
    [~, meta] = readKcatSource(fullfile(project, 'data', 'kcats', 'brenda.tsv'));
    verifyEqual(testCase, meta.parameters, '{"aggregation":"median","force_wildcard_level":0}')
    report = updateKcatSources(base, 'modelAdapter', adapter);
    verifyEqual(testCase, report.status', {'skipped','current','current'})
    % gecko_classic asks for the BRENDA maximum: brenda.tsv is outdated.
    adapter.params.kcat_assignment = struct('preset', 'gecko_classic');
    report = updateKcatSources(base, 'dryRun', true, 'modelAdapter', adapter);
    verifyEqual(testCase, report.status{2}, 'outdated')
end

function [base, adapter, project] = buildProject()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    model = getGeckoTestModel();
    base = makeEcModel(model, false, adapter);
    base = getECfromGEM(base);
    % DLKcat output for this model.
    rows = {'R2_EXP_1','G1','m1','C','MRAL','12.5'; 'R2_EXP_1','G2','m1','C','MNTD','6.5'; ...
            'R3','G4','m1','C','MDFM','3.25'};
    fid = fopen(fullfile(project, 'data', 'DLKcat.tsv'), 'w');
    fprintf(fid, 'rxns\tgenes\tsubstrates\tsmiles\tsequence\tkcat\n');
    for i = 1:size(rows, 1)
        fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', rows{i,:});
    end
    fclose(fid);
end

function writeLines(file, lines)
folder = fileparts(file);
if ~isfolder(folder)
    mkdir(folder);
end
fid = fopen(file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
