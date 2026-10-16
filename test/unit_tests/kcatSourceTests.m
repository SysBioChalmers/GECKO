%run this test case with the command
%results = runtests('kcatSourceTests.m')

function tests = kcatSourceTests
    tests = functiontests(localfunctions);
end

function testGetEcRxnKeys_tc0070(testCase)
    % Full and light models give the same conventional reactions and
    % enzyme sets.
    [full, light] = buildKcatModels();
    [baseF, enzF] = getEcRxnKeys(full);
    [baseL, enzL] = getEcRxnKeys(light);
    verifyEqual(testCase, baseF(strcmp(full.ec.rxns, 'R2_EXP_1')), {'R2'})
    verifyEqual(testCase, baseF(strcmp(full.ec.rxns, 'R2_REV_EXP_2')), {'R2_REV'})
    verifyEqual(testCase, sort(unique(strcat(baseF, '|', enzF))), sort(unique(strcat(baseL, '|', enzL))))
    for i = 1:numel(enzF)
        parts = strsplit(enzF{i}, ' + ');
        verifyEqual(testCase, parts, sort(parts))
    end
end

function testSaveLoadRoundTrip_tc0071(testCase)
    % What saveKcatSource writes, loadKcatSource maps back to the same
    % model.ec rows and values.
    [full, ~, adapter] = buildKcatModels();
    list = fuzzyKcatMatching(full, [], adapter);
    folder = tempname;
    file = saveKcatSource(full, list, 'brenda', 'folder', folder, ...
        'parameters', struct('aggregation','max'), 'data', 'BRENDA release test');
    [records, meta] = readKcatSource(file);
    verifyEqual(testCase, meta.source, 'brenda')
    verifyEqual(testCase, meta.kind, 'gathered')
    verifyEqual(testCase, meta.data, 'BRENDA release test')
    verifyEqual(testCase, meta.parameters, '{"aggregation":"max"}')
    keep = list.kcats > 0;
    verifyEqual(testCase, numel(records.kcat), sum(keep))
    loaded = loadKcatSource(full, file);
    verifyEqual(testCase, sort(loaded.rxns), sort(list.rxns(keep)))
    [~, a] = sort(list.rxns(keep)); [~, b] = sort(loaded.rxns);
    kept = find(keep);
    verifyEqual(testCase, loaded.kcats(b), list.kcats(kept(a)))
    verifyEqual(testCase, loaded.wildcardLvl(b), list.wildcardLvl(kept(a)))
    verifyEqual(testCase, loaded.origin(b), list.origin(kept(a)))
    verifyEqual(testCase, loaded.source, 'brenda')
    subs = list.substrates(kept(a));
    verifyEqual(testCase, loaded.substrates(b), cellfun(@(c) strjoin(c, '; '), subs, 'UniformOutput', false))
    verifyEqual(testCase, loaded.eccodes(b), list.eccodes(kept(a)))
    % Rows are sorted, so a rebuilt model writes the same file.
    sorted = sortrows([records.rxn records.enzymes records.substrate records.source]);
    verifyEqual(testCase, [records.rxn records.enzymes records.substrate records.source], sorted)
end

function testFullFileOnLightModel_tc0072(testCase)
    % A file written from the full model applies to the light model.
    [full, light, adapter] = buildKcatModels();
    list = fuzzyKcatMatching(full, [], adapter);
    file = saveKcatSource(full, list, 'brenda', 'folder', tempname);
    [loaded, report] = loadKcatSource(light, file);
    verifyEmpty(testCase, report.unmatched)
    [baseL, enzL] = getEcRxnKeys(light);
    [records] = readKcatSource(file);
    for r = 1:numel(records.kcat)
        hit = strcmp(baseL, records.rxn{r}) & strcmp(enzL, records.enzymes{r});
        verifyEqual(testCase, sort(loaded.rxns(loaded.line == records.line(r))), sort(light.ec.rxns(hit)))
    end
end

function testMatchingRules_tc0073(testCase)
    % Rows with only enzymes or only a reaction, unmatched and partial rows.
    [full, ~, adapter] = buildKcatModels();
    [base, enz] = getEcRxnKeys(full);
    multi = find(contains(enz, ' + '), 1);
    parts = strsplit(enz{multi}, ' + ');
    folder = tempname; mkdir(folder);
    file = fullfile(folder, 'custom.tsv');
    fid = fopen(file, 'w');
    fprintf(fid, '# source: custom\n# kind: curated\n');
    fprintf(fid, 'rxn\tenzymes\tsubstrate\tkcat\tlog10_sd\tsource\tdetail\tec_code\tnote\n');
    fprintf(fid, '\t%s\t\t10\t\tcustom\t\t\tenzymes only, reversed order\n', strjoin(fliplr(parts), ' + '));
    fprintf(fid, '%s\t\t\t20\t\tcustom\t\t\treaction only\n', base{multi});
    fprintf(fid, 'R999\tX1\t\t30\t\tcustom\t\t\tunmatched\n');
    fprintf(fid, '%s\t%s + XNEW\t\t40\t\tcustom\t\t\tpartial\n', base{multi}, enz{multi});
    fclose(fid);
    [loaded, report] = loadKcatSource(full, file);
    verifyEqual(testCase, sort(loaded.rxns(loaded.line == 4)), sort(full.ec.rxns(strcmp(enz, enz{multi}))))
    verifyEqual(testCase, sort(loaded.rxns(loaded.line == 5)), sort(full.ec.rxns(strcmp(base, base{multi}))))
    verifyEqual(testCase, report.unmatched.line, [6; 7])
    verifyEqual(testCase, unique(report.partial.line), 7)
    % Curated files are never written: by name, by 'kind: curated' and by
    % 'source: custom', in any case, and the file is left as it was.
    list = fuzzyKcatMatching(full, [], adapter);
    before = fileread(file);
    verifyError(testCase, @() saveKcatSource(full, list, 'custom', 'folder', folder), 'saveKcatSource:curated')
    verifyEqual(testCase, fileread(file), before)
    writeLines(fullfile(folder, 'mine.tsv'), {'# source: Custom', 'rxn	kcat	source'});
    verifyError(testCase, @() saveKcatSource(full, list, 'mine', 'folder', folder), 'saveKcatSource:curated')
    writeLines(fullfile(folder, 'mine2.tsv'), {'# kind: Curated', 'rxn	kcat	source'});
    verifyError(testCase, @() saveKcatSource(full, list, 'mine2', 'folder', folder), 'saveKcatSource:curated')
end

function testKcatsFromBrendaAndDLKcat_tc0074(testCase)
    % The source functions write their files into data/kcats.
    [full, ~, adapter] = buildKcatModels();
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    list = kcatsFromBrenda(full, 'aggregation', 'median', 'modelAdapter', adapter);
    [records, meta] = readKcatSource(fullfile(project, 'data', 'kcats', 'brenda.tsv'));
    verifyEqual(testCase, meta.generator(1:15), 'kcatsFromBrenda')
    verifyEqual(testCase, meta.parameters, '{"aggregation":"median","force_wildcard_level":0}')
    verifyEqual(testCase, meta.data, 'BRENDA release test-fixture')
    verifyEqual(testCase, numel(records.kcat), sum(list.kcats > 0))

    % DLKcat output for this model, read without running DLKcat.
    rows = {'R2_EXP_1','G1','m1','C','MRAL','12.5'; 'R3','G4','m1','C','MDFM','3.25'};
    fid = fopen(fullfile(project, 'data', 'DLKcat.tsv'), 'w');
    fprintf(fid, 'rxns\tgenes\tsubstrates\tsmiles\tsequence\tkcat\n');
    for i = 1:size(rows, 1)
        fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', rows{i,:});
    end
    fclose(fid);
    kcatsFromDLKcat(full, 'modelAdapter', adapter);
    [records, meta] = readKcatSource(fullfile(project, 'data', 'kcats', 'dlkcat.tsv'));
    verifyEqual(testCase, meta.source, 'dlkcat')
    verifyEqual(testCase, sort(records.kcat), [3.25; 12.5])
    verifyEqual(testCase, unique(records.source), {'DLKcat'})
end

function testKcatsFromOKPStored_tc0075(testCase)
    % Stored OKP results: a method's own result of its current job, or a
    % result stored by fetchOpenKineticsPredictor alone when its predictions
    % are of the requested method.
    [full, ~, adapter] = buildKcatModels();
    full.metSmiles = repmat({''}, numel(full.mets), 1);
    full.metSmiles(strcmp(full.mets, 'm1c')) = {'C(C1C)O'};
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    data = fullfile(project, 'data');
    okpCsv = {'kcat (1/s),Source kcat,Extra Info kcat,Protein Sequence,Substrate', ...
              '12.5,Prediction from CataPro,,MRAL,C(C1C)O'};
    writeLines(fullfile(data, 'OKP_output.csv'), okpCsv);
    writeLines(fullfile(data, 'OKP_job.txt'), {'jobId: abc123', 'method: CataPro'});
    [list, done] = kcatsFromOKP(full, 'method', "CataPro", 'useStored', true, 'modelAdapter', adapter);
    verifyTrue(testCase, done)
    [records, meta] = readKcatSource(fullfile(data, 'kcats', 'okp_catapro.tsv'));
    verifyEqual(testCase, meta.data, 'OKP job abc123')
    verifyEqual(testCase, numel(records.kcat), sum(list.kcats > 0))
    verifyEqual(testCase, unique(records.source), {'OKP-CataPro'})
    verifyTrue(testCase, isfile(fullfile(data, 'OKP_output_catapro.csv')))

    % The shared job file names UniKP, but the shared result holds CataPro
    % predictions (fetched after UniKP was submitted): refused.
    writeLines(fullfile(data, 'OKP_job.txt'), {'jobId: uni1', 'method: UniKP'});
    verifyError(testCase, @() kcatsFromOKP(full, 'method', 'UniKP', 'useStored', true, ...
        'modelAdapter', adapter), 'kcatsFromOKP:noStored')
    verifyFalse(testCase, isfile(fullfile(data, 'kcats', 'okp_unikp.tsv')))
    % A method's stored result of an earlier job is refused too.
    writeLines(fullfile(data, 'OKP_job_kinform_h.txt'), {'jobId: new2', 'method: KinForm-H', 'resultOf: old1'});
    writeLines(fullfile(data, 'OKP_output_kinform_h.csv'), okpCsv);
    verifyError(testCase, @() kcatsFromOKP(full, 'method', 'KinForm-H', 'useStored', true, ...
        'modelAdapter', adapter), 'kcatsFromOKP:noStored')
    % ecRxns when no job is submitted, and a misspelled option, are errors.
    verifyError(testCase, @() kcatsFromOKP(full, 'method', 'CataPro', 'useStored', true, ...
        'ecRxns', true(numel(full.ec.rxns), 1), 'modelAdapter', adapter), 'kcatsFromOKP:ecRxns')
    verifyError(testCase, @() kcatsFromOKP(full, 'method', 'UniKP', 'wiat', true, ...
        'modelAdapter', adapter), 'kcatsFromOKP:arguments')
end

function testGuardsNumbersAndJson_tc0078(testCase)
    [full, ~, adapter] = buildKcatModels();
    project = tempname;
    copyfile(fullfile(adapter.params.path, 'data'), fullfile(project, 'data'));
    adapter.params.path = project;
    folder = fullfile(project, 'data', 'kcats');
    mkdir(folder);
    % The curated guard reads quoted header lines and ignores the case of keys.
    writeLines(fullfile(folder, 'a.tsv'), {'"# kind: curated"', 'rxn	kcat	source'});
    writeLines(fullfile(folder, 'b.tsv'), {'# Kind: Curated', 'rxn	kcat	source'});
    writeLines(fullfile(folder, 'c.tsv'), {'# SOURCE: custom', 'rxn	kcat	source'});
    list = fuzzyKcatMatching(full, [], adapter);
    for f = {'a','b','c'}
        verifyError(testCase, @() saveKcatSource(full, list, f{1}, 'modelAdapter', adapter), 'saveKcatSource:curated')
    end
    % The source functions check the file before they gather.
    verifyError(testCase, @() kcatsFromBrenda(full, 'saveAs', "custom", 'modelAdapter', adapter), 'saveKcatSource:curated')
    verifyError(testCase, @() kcatsFromDLKcat(full, 'ecRxns', true(numel(full.ec.rxns), 1), ...
        'modelAdapter', adapter), 'kcatsFromDLKcat:ecRxns')

    % Infinite, NaN and non-positive kcats are not written, and not applied.
    list.kcats(1) = Inf; list.kcats(2) = NaN; list.kcats(3) = -2;
    file = saveKcatSource(full, list, "numbers", 'modelAdapter', adapter, 'parameters', ...
        struct('z', 1, 'a', struct('y', ['caf' char(233)], 'b', true), 'n', NaN, ...
        'f', 1/3, 'g', 1e-7, 'h', 1e15, 'k', 1e20, 'l', {{5}}, 's', 'tab	quote"'));
    [records, meta] = readKcatSource(file);
    verifyTrue(testCase, all(isfinite(records.kcat) & records.kcat > 0))
    % The same text as geckopy's test_parameters_json_matches_gecko.
    verifyEqual(testCase, meta.parameters, ['{"a":{"b":true,"y":"caf' char(233) '"},"f":0.3333333333333333,' ...
        '"g":1e-07,"h":1e+15,"k":1e+20,"l":[5],"n":null,"s":"tab\tquote\"","z":1}'])
    writeLines(fullfile(folder, 'inf.tsv'), {'rxn	enzymes	kcat	source', 'R3	P4	Inf	custom', 'R3	P4	1e999	custom'});
    [loaded, report] = loadKcatSource(full, fullfile(folder, 'inf.tsv'));
    verifyEmpty(testCase, loaded.rxns)
    verifyEqual(testCase, report.invalid.line, [2; 3])

    % Rows are sorted by every column, in the same order as geckopy.
    writeLines(fullfile(folder, 'unsorted.tsv'), {'rxn	enzymes	substrate	kcat	source	detail', ...
        'R3	P4	m1	2	b	', 'R3	P4	m1	2	a	wc=1', 'R3	P4		5	a	', 'R2_REV	P3	m2	1	a	'});
    l2 = loadKcatSource(full, fullfile(folder, 'unsorted.tsv'));
    file = saveKcatSource(full, l2, 'sorted', 'modelAdapter', adapter);
    records = readKcatSource(file);
    verifyEqual(testCase, records.rxn', {'R2_REV','R3','R3','R3'})
    verifyEqual(testCase, records.substrate', {'m2','','m1','m1'})
    verifyEqual(testCase, records.source', {'a','a','a','b'})
end

function testEnzymelessRowsAndLightToFull_tc0076(testCase)
    % Rows of reactions without enzymes are not written, so they cannot
    % spread to their isozymes; a light-model file applies to the full model.
    [full, light, adapter] = buildKcatModels();
    list = fuzzyKcatMatching(full, [], adapter);
    i = find(strcmp(full.ec.rxns, 'R2_EXP_2'));
    full.ec.rxnEnzMat(i, :) = 0;
    file = saveKcatSource(full, list, 'brenda', 'folder', tempname);
    loaded = loadKcatSource(full, file);
    verifyFalse(testCase, any(strcmp(loaded.rxns, 'R2_EXP_2')))
    verifyEqual(testCase, sum(strcmp(loaded.rxns, 'R2_EXP_1')), 1)

    listL = fuzzyKcatMatching(light, [], adapter);
    fileL = saveKcatSource(light, listL, 'brenda', 'folder', tempname);
    fullAgain = buildKcatModels();
    [~, report] = loadKcatSource(fullAgain, fileL);
    verifyEmpty(testCase, report.unmatched)
end

function testSharedFixture_tc0077(testCase)
    % The fixture file, read the same way by geckopy: byte-order mark, CRLF,
    % a quoted header line and field, comments and blank lines between
    % rows, a row of tabs, a decimal comma, a negative kcat.
    file = fullfile(findGECKOroot, 'test', 'unit_tests', 'ecTestGEM', 'data', 'kcatSourceFixture.tsv');
    [records, meta] = readKcatSource(file);
    verifyEqual(testCase, meta.source, 'fixture')
    verifyEqual(testCase, meta.parameters, '{"aggregation":"median"}')
    verifyEqual(testCase, meta.kind, 'curated')
    verifyEqual(testCase, records.line', [7 9 11 12 13 14 15])
    verifyEqual(testCase, records.rxn', {'R2','R3','R5','R2_REV','','R9','R9'})
    verifyEqual(testCase, records.substrate{2}, '2,3-bisphospho-D-glycerate')
    verifyEqual(testCase, records.kcat', [12.5 1e-05 NaN 7 3 -1 4])
    verifyEqual(testCase, records.log10_sd(1), 0.3)

    [full] = buildKcatModels();
    [loaded, report] = loadKcatSource(full, file);
    verifyEqual(testCase, loaded.rxns(loaded.line == 7), {'R2_EXP_1'})
    verifyEqual(testCase, sort(loaded.rxns(loaded.line == 12)), {'R2_REV_EXP_1';'R2_REV_EXP_2'})
    verifyEqual(testCase, sort(loaded.rxns(loaded.line == 13)), {'R2_EXP_2';'R2_REV_EXP_2'})
    verifyEqual(testCase, report.invalid.line, [11; 14])
    verifyEqual(testCase, report.unmatched.line, 15)

    % A row with more filled fields than columns is an error.
    bad = [tempname '.tsv'];
    writeLines(bad, {'rxn	kcat	source', 'R2	1	brenda	extra'});
    verifyError(testCase, @() readKcatSource(bad), 'readKcatSource:fields')
end

function writeLines(file, lines)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end

function [full, light, adapter] = buildKcatModels()
    geckoPath = findGECKOroot;
    adapter = ModelAdapterManager.getAdapter(fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'));
    model = getGeckoTestModel();
    full = makeEcModel(model, false, adapter);
    full = getECfromGEM(full);
    light = makeEcModel(model, true, adapter);
    light = getECfromGEM(light);
end
