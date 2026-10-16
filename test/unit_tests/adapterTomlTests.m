%run this test case with the command
%results = runtests('adapterTomlTests.m')

function tests = adapterTomlTests
    tests = functiontests(localfunctions);
end

function testReadAdapterToml_tc0061(testCase)
    % Every TOML feature that model adapter files use.
    file = writeTemp({ ...
        '# comment', ...
        'name = "a \"quoted\" # not a comment"  # comment', ...
        'literal = ''C:\path''', ...
        'int = 559_292', ...
        'float = 2.5e-1', ...
        'yes = true', ...
        'list = ["a", "b",', ...
        '        "c"]  # over two lines', ...
        'nums = [0.4, 0.2, 0.1]', ...
        'empty = []', ...
        '', ...
        '[kegg]', ...
        'id = "sce"', ...
        '', ...
        '[evotune.source_groups.okp]', ...
        'sources = ["OKP"]', ...
        'match_okp = true', ...
        '', ...
        '[[stage]]', ...
        'action = "select"', ...
        '[[stage]]', ...
        'action = "fill"', ...
        'sources = ["isozymes"]'});
    d = readAdapterToml(file);
    verifyEqual(testCase, d.name, 'a "quoted" # not a comment')
    verifyEqual(testCase, d.literal, 'C:\path')
    verifyEqual(testCase, d.int, 559292)
    verifyEqual(testCase, d.float, 0.25)
    verifyEqual(testCase, d.yes, true)
    verifyEqual(testCase, d.list, {'a','b','c'})
    verifyEqual(testCase, d.nums, [0.4 0.2 0.1])
    verifyEqual(testCase, d.empty, {})
    verifyEqual(testCase, d.kegg.id, 'sce')
    verifyEqual(testCase, d.evotune.source_groups.okp.sources, {'OKP'})
    verifyTrue(testCase, d.evotune.source_groups.okp.match_okp)
    verifyEqual(testCase, numel(d.stage), 2)
    verifyEqual(testCase, d.stage(2).action, 'fill')
    verifyEqual(testCase, d.stage(2).sources, {'isozymes'})
end

function testReadAdapterTomlErrors_tc0062(testCase)
    % Unsupported TOML stops with the line number.
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = { a = 1 }'})), 'GECKO:readAdapterToml')
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = 2026-10-15'})), 'GECKO:readAdapterToml')
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = unquoted'})), 'GECKO:readAdapterToml')
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = 1', 'x = 2'})), 'GECKO:readAdapterToml')
    try
        readAdapterToml(writeTemp({'a = 1', 'b = """x"""'}));
    catch err
        verifySubstring(testCase, err.message, 'line 2')
    end
end

function testWriteAdapterToml_tc0063(testCase)
    % What writeAdapterToml writes, readAdapterToml reads back unchanged.
    s.text = 'with "quotes" and \ backslash';
    s.number = 0.1;
    s.integer = 150;
    s.flag = false;
    s.list = {'dlkcat','brenda'};
    s.vector = [0.4 0.2];
    s.table.inner.value = 1e-8;
    file = [tempname '.toml'];
    writeAdapterToml(s, file, 'header');
    verifyEqual(testCase, readAdapterToml(file), s)
end

function testLoadAdapterParamsDefaults_tc0064(testCase)
    % A minimal file gets geckopy's defaults, under GECKO's field names.
    folder = tempname; mkdir(folder);
    writeTemp({'conv_gem = "models/x.xml"', 'org_name = "testus testus"', 'p_tot = 0.4'}, ...
        fullfile(folder, 'model_adapter.toml'));
    params = loadAdapterParams(folder);
    verifyEqual(testCase, params.path, folder)
    verifyEqual(testCase, params.convGEM, fullfile(folder, 'models/x.xml'))
    verifyEqual(testCase, params.Ptot, 0.4)
    verifyEqual(testCase, params.sigma, 0.5)
    verifyEqual(testCase, params.kegg.ID, 'sce')
    verifyEqual(testCase, params.kcatAggregation, 'max')
    verifyEqual(testCase, params.evotune.kcatSources, {'dlkcat','brenda','custom'})
    verifyEqual(testCase, params.evotune.sigma0logSource, [0.4; 0.2; 0.1])
    verifyEqual(testCase, params.evotune.maxGenerations, 150)
    verifyEqual(testCase, params.okp.handleLongSequences, 'truncate')
    verifyEmpty(testCase, params.complex.taxonomicID)
    % A misspelled key is an error, as in geckopy.
    writeTemp({'conv_gem = "x.xml"', 'org_name = "x"', 'sigmaa = 1'}, fullfile(folder, 'model_adapter.toml'));
    verifyError(testCase, @() loadAdapterParams(folder), 'GECKO:loadAdapterParams')
end

function testLoadAdapterParamsSourceGroups_tc0065(testCase)
    % Trust tiers with several sources and the OKP method.
    folder = tempname; mkdir(folder);
    writeTemp({'conv_gem = "x.xml"', 'org_name = "x"', ...
        '[okp]', 'method = "EITLEM"', ...
        '[evotune.source_groups.okp]', 'sources = ["OpenKineticsPredictor"]', 'match_okp = true', ...
        '[evotune.source_groups.database]', 'sources = ["brenda", "okp_brenda"]', ...
        '[evotune.sigma0_log_source]', 'okp = 0.25', 'database = 0.2', ...
        '[matlab]', 'make_model_anaerobic = "anaerobicModel_GECKO"', ...
        '[matlab.params]', 'c_UptakeExp = 0.64'}, fullfile(folder, 'model_adapter.toml'));
    params = loadAdapterParams(fullfile(folder, 'model_adapter.toml'));
    verifyEqual(testCase, params.evotune.kcatSources, {'OpenKineticsPredictor','EITLEM','brenda','okp_brenda'})
    verifyEqual(testCase, params.evotune.sigma0logSource, [0.25; 0.25; 0.2; 0.2])
    verifyEqual(testCase, params.matlab.makeModelAnaerobic, 'anaerobicModel_GECKO')
    verifyEqual(testCase, params.c_UptakeExp, 0.64)
end

function testConvertAdapterToToml_tc0066(testCase)
    % Converting an adapter class and loading the TOML gives the same
    % parameters back.
    geckoPath = findGECKOroot;
    legacyFolder = tempname; mkdir(legacyFolder);
    writeTemp(legacyAdapterClass(), fullfile(legacyFolder, 'LegacyFixtureAdapter.m'));
    adapters = {fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'), ...
                fullfile(legacyFolder,'LegacyFixtureAdapter.m')};
    for k = 1:numel(adapters)
        original = ModelAdapterManager.getAdapter(adapters{k});
        folder = tempname; mkdir(folder);
        file = fullfile(folder, 'model_adapter.toml');
        report = convertAdapterToToml(original, file);
        verifyTrue(testCase, any(strcmp(report.methods, 'getSpontaneousReactions')))
        converted = TomlAdapter(file);
        expected = original.params;
        actual = converted.params;
        % The project folder is where the TOML file is; the model path is
        % relative to it.
        [~, convName] = fileparts(expected.convGEM);
        [~, convName2] = fileparts(actual.convGEM);
        verifyEqual(testCase, convName2, convName)
        expected = rmfield(expected, {'path','convGEM'});
        verifyEqual(testCase, keepFields(actual, expected), expected)
    end
    verifyEqual(testCase, report.unmapped, {'c_UptakeExp'})
end

function lines = legacyAdapterClass()
% An adapter class as written for GECKO 3 and 4: evotune trust tiers and a
% parameter that geckopy does not know.
lines = { ...
    'classdef LegacyFixtureAdapter < ModelAdapter', ...
    '    methods', ...
    '        function obj = LegacyFixtureAdapter()', ...
    '            obj.params.path = fileparts(mfilename(''fullpath''));', ...
    '            obj.params.convGEM = fullfile(obj.params.path,''models'',''human-GEM.yml'');', ...
    '            obj.params.sigma = 0.1;', ...
    '            obj.params.Ptot = 0.505717472;', ...
    '            obj.params.f = 0.412;', ...
    '            obj.params.gR_exp = 0.020663429;', ...
    '            obj.params.c_UptakeExp = 0.641339301;', ...
    '            obj.params.org_name = ''homo sapiens'';', ...
    '            obj.params.kcatAggregation = ''median'';', ...
    '            obj.params.complex.taxonomicID = 9606;', ...
    '            obj.params.kegg.ID = ''hsa'';', ...
    '            obj.params.kegg.geneID = ''Ensembl'';', ...
    '            obj.params.uniprot.type = ''proteome'';', ...
    '            obj.params.uniprot.ID = ''UP000005640'';', ...
    '            obj.params.uniprot.geneIDfield = ''gene_primary'';', ...
    '            obj.params.uniprot.reviewed = true;', ...
    '            obj.params.c_source = ''MAR09034'';', ...
    '            obj.params.bioRxn = ''MAR13082'';', ...
    '            obj.params.enzyme_comp = ''Cytosol'';', ...
    '            obj.params.evotune.sigma0logDefault = 0.5;', ...
    '            obj.params.evotune.kcatSources = {''dlkcat'',''brenda'',''custom''};', ...
    '            obj.params.evotune.sigma0logSource = [0.4; 0.2; 0.1];', ...
    '            obj.params.evotune.maxGrowthWeight = 2.0;', ...
    '            obj.params.evotune.priorPenaltyWeight = 0.03;', ...
    '            obj.params.evotune.tieIsozymes = true;', ...
    '            obj.params.evotune.rmseThreshold = 0.2;', ...
    '            obj.params.evotune.maxGenerations = 212;', ...
    '        end', ...
    '        function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)', ...
    '            spont = contains(model.rxnNames,''spontaneous'');', ...
    '            spontRxnNames = model.rxnNames(spont);', ...
    '        end', ...
    '    end', ...
    'end'};
end

function testTomlAdapterSpontaneous_tc0067(testCase)
    folder = tempname; mkdir(fullfile(folder,'data'));
    writeTemp({'rxns	spontaneous', 'R1	0', 'R2	1', 'R3	1'}, fullfile(folder,'data','spont.tsv'));
    writeTemp({'conv_gem = "x.xml"', 'org_name = "x"', ...
        '[spontaneous]', 'name_contains = "spontaneous"', 'file = "data/spont.tsv"'}, ...
        fullfile(folder, 'model_adapter.toml'));
    adapter = ModelAdapterManager.getAdapter(folder);
    verifyClass(testCase, adapter, 'TomlAdapter')
    model.rxns = {'R1';'R2';'R4';'R5'};
    model.rxnNames = {'a';'b';'spontaneous c';'d'};
    [spont, names] = adapter.getSpontaneousReactions(model);
    verifyEqual(testCase, spont, [false; true; true; false])
    verifyEqual(testCase, names, {'R2';'R4'})
end

function testCreateGECKOprojectToml_tc0068(testCase)
    % A new project gets model_adapter.toml, which loads as it is.
    here = pwd;
    restore = onCleanup(@() cd(here));
    root = tempname; mkdir(root);
    createGECKOproject('name', 'ecTestProject', 'path', root);
    folder = fullfile(root, 'ecTestProject');
    verifyTrue(testCase, isfile(fullfile(folder, 'model_adapter.toml')))
    adapter = ModelAdapterManager.getAdapter(folder);
    verifyEqual(testCase, adapter.params.org_name, 'genus species')
    verifyEqual(testCase, adapter.params.evotune.kcatSources, {'dlkcat','brenda','custom'})
    verifyEqual(testCase, adapter.params.spontaneous.name_contains, 'spontaneous')
end

function s = keepFields(s, template)
% s reduced to the fields (recursively) that template has.
fields = fieldnames(s);
for k = 1:numel(fields)
    if ~isfield(template, fields{k})
        s = rmfield(s, fields{k});
    elseif isstruct(s.(fields{k})) && isstruct(template.(fields{k}))
        s.(fields{k}) = keepFields(s.(fields{k}), template.(fields{k}));
    end
end
s = orderfields(s, template);
end

function file = writeTemp(lines, file)
if nargin < 2
    file = [tempname '.toml'];
end
fid = fopen(file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
