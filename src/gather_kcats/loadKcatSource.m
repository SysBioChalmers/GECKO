function [kcatList, report] = loadKcatSource(model, name, varargin)
% loadKcatSource  Read a kcat source file and map its rows to an ecModel.
%
% Reads data/kcats/<name>.tsv (see saveKcatSource) and maps every row to
% the model.ec reactions it applies to, by conventional reaction and
% enzyme set (see getEcRxnKeys):
%
% 1. A row with rxn and enzymes applies to the model.ec rows with that
%    conventional reaction and that enzyme set.
% 2. A row with only enzymes applies to every model.ec row with that
%    enzyme set; a row with only rxn to every model.ec row of that
%    conventional reaction (curated rows; saveKcatSource always writes
%    both).
% 3. A row whose kcat is missing, not a number, infinite, or not above 0
%    is reported as invalid and not applied.
% 4. A row that matches nothing is reported, never dropped silently. If a
%    model.ec row has at least half of the row's enzymes, that partial
%    match is reported as well, but not applied.
%
% The enzymes of a row may be listed in any order.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% name : char
%     name of the source (data/kcats/<name>.tsv in the project folder), or
%     the path to a kcat source file.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     kcat values per model.ec reaction, as used by mergeKcats and
%     assignKcatValues: source (the source name), rxns, kcats,
%     substrates, eccodes, wildcardLvl and origin (from the detail
%     column), kcatSource (the per-row source), log10sd, and line (the
%     line of the file each entry comes from).
% report : struct
%     - invalid: table of the file rows with a missing or invalid kcat;
%     - unmatched: table of the file rows that apply to no model.ec row;
%     - partial: table of rows that only partly match an enzyme set, with
%       the model.ec reaction and its enzymes.
%
% Examples
% --------
%     kcatList = loadKcatSource(ecModel, 'brenda');
%     [kcatList, report] = loadKcatSource(ecModel, 'custom');
%
% See also
% --------
% saveKcatSource, readKcatSource, kcatsFromBrenda, getEcRxnKeys

checkNameValue(varargin, {'modelAdapter'}, 'loadKcatSource');
p = parseGECKOargs(varargin, {'modelAdapter', []});

if isfile(name)
    filename = name;
else
    modelAdapter = p.modelAdapter;
    if isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    filename = fullfile(modelAdapter.params.path, 'data', 'kcats', [name '.tsv']);
    if ~isfile(filename)
        error('loadKcatSource:missing', ['No kcat source file %s. Make it with the ' ...
            'kcatsFrom<Source> function of this source (kcatsFromBrenda, kcatsFromDLKcat, ' ...
            'kcatsFromOKP), or with saveKcatSource.'], filename)
    end
end
[records, meta] = readKcatSource(filename);
sourceName = '';
if isfield(meta, 'source')
    sourceName = meta.source;
end

[baseRxns, enzymeSets] = getEcRxnKeys(model);
nEc = numel(baseRxns);
nRec = numel(records.kcat);

% Lookup tables: key -> model.ec rows.
byPair = groupIndex(strcat(baseRxns, '|', enzymeSets));
byBase = groupIndex(baseRxns);
byEnz  = groupIndex(enzymeSets);
% accession -> model.ec rows containing it, for partial matches.
accRows = cell(nEc, 1);
for i = 1:nEc
    if ~isempty(enzymeSets{i})
        accRows{i} = strsplit(enzymeSets{i}, ' + ');
    end
end
accIdx = repelem((1:nEc)', cellfun(@numel, accRows));
byAcc = groupIndex([accRows{:}]', accIdx);

invalid = ~(records.kcat > 0) | ~isfinite(records.kcat);
matchRec = cell(nRec, 1);
matchEc  = cell(nRec, 1);
unmatched = false(nRec, 1);
partialRec = zeros(0, 1); partialEc = zeros(0, 1);
for r = find(~invalid)'
    rxn = records.rxn{r};
    enz = normalizeEnzymes(records.enzymes{r});
    if ~isempty(rxn) && ~isempty(enz)
        hit = lookup(byPair, [rxn '|' enz]);
    elseif ~isempty(enz)
        hit = lookup(byEnz, enz);
    elseif ~isempty(rxn)
        hit = lookup(byBase, rxn);
    else
        hit = [];
    end
    if isempty(hit)
        unmatched(r) = true;
        if ~isempty(enz)
            wanted = strsplit(enz, ' + ');
            candidates = unique(cell2mat(cellfun(@(a) lookup(byAcc, a), wanted(:), 'UniformOutput', false)));
            if ~isempty(rxn)
                candidates = intersect(candidates, lookup(byBase, rxn));
            end
            for c = candidates(:)'
                if mean(ismember(wanted, accRows{c})) >= 0.5
                    partialRec(end+1, 1) = r; %#ok<AGROW>
                    partialEc(end+1, 1) = c; %#ok<AGROW>
                end
            end
        end
    else
        matchRec{r} = repmat(r, numel(hit), 1);
        matchEc{r} = hit(:);
    end
end
rec = vertcat(matchRec{:}, zeros(0, 1));
ec  = vertcat(matchEc{:}, zeros(0, 1));

kcatList.source      = sourceName;
kcatList.rxns        = reshape(model.ec.rxns(ec), [], 1);
kcatList.kcats       = records.kcat(rec);
kcatList.substrates  = records.substrate(rec);
kcatList.eccodes     = records.ec_code(rec);
kcatList.kcatSource  = records.source(rec);
kcatList.log10sd     = records.log10_sd(rec);
[kcatList.wildcardLvl, kcatList.origin] = parseDetail(records.detail(rec));
kcatList.line        = records.line(rec);

report.invalid = table(records.line(invalid), records.rxn(invalid), records.enzymes(invalid), ...
    records.kcat(invalid), 'VariableNames', {'line','rxn','enzymes','kcat'});
report.unmatched = table(records.line(unmatched), records.rxn(unmatched), records.enzymes(unmatched), ...
    records.kcat(unmatched), 'VariableNames', {'line','rxn','enzymes','kcat'});
report.partial = table(records.line(partialRec), records.rxn(partialRec), records.enzymes(partialRec), ...
    reshape(model.ec.rxns(partialEc), [], 1), enzymeSets(partialEc), ...
    'VariableNames', {'line','rxn','enzymes','ecRxn','ecEnzymes'});
if any(invalid)
    printOrange(sprintf('%d of %d rows of %s have no valid kcat; see report.invalid.\n', ...
        sum(invalid), nRec, filename));
end
if any(unmatched)
    printOrange(sprintf('%d of %d rows of %s apply to no reaction of this model; see report.unmatched.\n', ...
        sum(unmatched), nRec, filename));
end
end

function map = groupIndex(keys, rows)
% containers.Map from each key to the rows (default 1:n) that have it.
if nargin < 2
    rows = (1:numel(keys))';
end
map = containers.Map('KeyType', 'char', 'ValueType', 'any');
if isempty(keys)
    return
end
[u, ~, g] = unique(keys(:));
groups = accumarray(g, rows(:), [], @(x) {sort(x)});
for i = 1:numel(u)
    if ~isempty(u{i})
        map(u{i}) = groups{i};
    end
end
end

function rows = lookup(map, key)
if isKey(map, key)
    rows = map(key);
else
    rows = zeros(0, 1);
end
end

function enz = normalizeEnzymes(enz)
% Sorted and joined as getEcRxnKeys does, so hand-written rows match too.
if isempty(enz)
    return
end
parts = strtrim(strsplit(enz, '+'));
parts = unique(parts(~cellfun(@isempty, parts)));
enz = strjoin(sort(parts(:)'), ' + ');
end

function [wc, origin] = parseDetail(detail)
wc = nan(numel(detail), 1);
origin = nan(numel(detail), 1);
for i = 1:numel(detail)
    t = regexp(detail{i}, 'wc=(\d+)', 'tokens', 'once');
    if ~isempty(t)
        wc(i) = str2double(t{1});
    end
    t = regexp(detail{i}, 'origin=(\d+)', 'tokens', 'once');
    if ~isempty(t)
        origin(i) = str2double(t{1});
    end
end
end
