function filename = saveKcatSource(model, kcatList, name, varargin)
% saveKcatSource  Write a kcatList to a kcat source file of the project.
%
% Writes the kcat values of kcatList to data/kcats/<name>.tsv in the
% project folder, the file that loadKcatSource reads back. Rows are keyed
% by their conventional reaction and enzyme set (see getEcRxnKeys), not
% by model.ec.rxns, so the file stays valid when the ecModel is rebuilt.
% kcatsFromBrenda, kcatsFromDLKcat and kcatsFromOKP write these files for
% their sources; use this function for a kcatList made in another way.
%
% The file is UTF-8 text. It starts with '#' header lines (source,
% generator, parameters, data, kind, created), then a tab-separated table
% with the columns rxn, enzymes, substrate, kcat, log10_sd, source,
% detail, ec_code and note, sorted by all columns from rxn onwards.
% Left out are:
% - rows whose kcat is not a positive finite number (0 means no value
%   was found);
% - rows of model.ec reactions without enzymes, which have no enzyme
%   constraint and could not be told apart from their isozymes.
%
% Parameters
% ----------
% model : struct
%     the ecModel the kcatList refers to.
% kcatList : struct
%     kcat values as returned by kcatsFromBrenda, fuzzyKcatMatching,
%     readDLKcatOutput, readOpenKineticsPredictorOutput or mergeKcats:
%     fields rxns and kcats, optionally source, substrates, eccodes,
%     wildcardLvl, origin, kcatSource and log10sd.
% name : char
%     name of the source; the file is <name>.tsv.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% folder : char
%     folder to write to (default: data/kcats in the project folder).
% generator : char
%     what produced the values (default: 'saveKcatSource').
% parameters : struct
%     settings used to produce the values, written as compact JSON with
%     sorted keys (default: none).
% data : char
%     version of the data the values come from, e.g. a BRENDA release
%     (default: empty).
% kind : char
%     'gathered' (default) or 'curated'.
%
% A file whose header says 'kind: curated' or 'source: custom', or named
% custom.tsv, is never written.
%
% Returns
% -------
% filename : char
%     the file written.
%
% Examples
% --------
%     kcatList = fuzzyKcatMatching(ecModel);
%     saveKcatSource(ecModel, kcatList, 'mybrenda');
%
% See also
% --------
% kcatsFromBrenda, loadKcatSource, readKcatSource, getEcRxnKeys

names = {'modelAdapter','folder','generator','parameters','data','kind'};
checkNameValue(varargin, names, 'saveKcatSource');
p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'folder',       []; ...
    'generator',    'saveKcatSource'; ...
    'parameters',   struct(); ...
    'data',         ''; ...
    'kind',         'gathered'});
kind = validatestring(p.kind, {'gathered','curated'});
name = char(name);
filename = kcatSourceFile(name, p.folder, p.modelAdapter);

[baseRxns, enzymeSets] = getEcRxnKeys(model);
[found, idx] = ismember(kcatList.rxns(:), model.ec.rxns);
if ~all(found)
    error('saveKcatSource:unknownRxns', '%s', ravenList( ...
        'These kcatList reactions are not in model.ec.rxns:', kcatList.rxns(~found), false));
end

n = numel(kcatList.rxns);
kcats = double(kcatList.kcats(:));
if numel(kcats) ~= n
    error('saveKcatSource:fields', 'kcatList.kcats and kcatList.rxns differ in length.')
end
substrate = textColumn(kcatList, 'substrates', n);
ecCode    = textColumn(kcatList, 'eccodes', n);
if isfield(kcatList, 'kcatSource')
    source = textColumn(kcatList, 'kcatSource', n);
elseif isfield(kcatList, 'source')
    source = repmat({char(kcatList.source)}, n, 1);
else
    source = repmat({name}, n, 1);
end
wc     = numColumn(kcatList, 'wildcardLvl', n);
origin = numColumn(kcatList, 'origin', n);
sd     = numColumn(kcatList, 'log10sd', n);

keep = kcats > 0 & isfinite(kcats);
noEnzymes = keep & cellfun(@isempty, enzymeSets(idx));
if any(noEnzymes)
    fprintf('%d kcats of reactions without enzymes are not written to %s.\n', sum(noEnzymes), filename);
end
keep = find(keep & ~noEnzymes);

clean = @(c) regexprep(string(c(:)), '[\t\r\n]', ' ');
rows = strings(numel(keep), 9);
if ~isempty(keep)
    wcText = strings(numel(keep), 1);
    hasWc = ~isnan(wc(keep));
    wcText(hasWc) = "wc=" + string(fix(wc(keep(hasWc))));
    originText = strings(numel(keep), 1);
    hasOrigin = ~isnan(origin(keep));
    originText(hasOrigin) = "origin=" + string(fix(origin(keep(hasOrigin))));
    detail = wcText + "; " + originText;
    detail(~hasWc) = originText(~hasWc);
    detail(~hasOrigin) = wcText(~hasOrigin);
    rows(:, 1) = clean(baseRxns(idx(keep)));
    rows(:, 2) = clean(enzymeSets(idx(keep)));
    rows(:, 3) = clean(substrate(keep));
    rows(:, 4) = formatNumbers(kcats(keep));
    rows(:, 5) = formatNumbers(sd(keep));
    rows(:, 6) = clean(source(keep));
    rows(:, 7) = detail;
    rows(:, 8) = clean(ecCode(keep));
    % Sorted by the columns in this order, as geckopy sorts tuples: one key
    % per row, the fields joined by a character below any in the text.
    [~, order] = sort(join(rows(:, [1 2 3 6 4 5 7 8 9]), char(1), 2));
    rows = rows(order, :);
end

folder = fileparts(filename);
if ~isfolder(folder)
    mkdir(folder);
end
fid = fopen(filename, 'w', 'n', 'UTF-8');
if fid < 0
    error('saveKcatSource:write', 'Cannot write %s', filename)
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# source: %s\n', name);
fprintf(fid, '# generator: %s (GECKO %s)\n', cleanText(p.generator), geckoVersion());
fprintf(fid, '# parameters: %s\n', kcatSourceJson(p.parameters));
fprintf(fid, '# data: %s\n', cleanText(p.data));
fprintf(fid, '# kind: %s\n', kind);
fprintf(fid, '# created: %s\n', char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
fprintf(fid, 'rxn\tenzymes\tsubstrate\tkcat\tlog10_sd\tsource\tdetail\tec_code\tnote\n');
if ~isempty(rows)
    fprintf(fid, '%s\n', join(rows, sprintf('\t'), 2));
end
end

function col = textColumn(s, field, n)
% A kcatList field as an n x 1 cell of char; a cell of several entries
% (e.g. all substrates of a reaction) is joined by '; '.
col = repmat({''}, n, 1);
if ~isfield(s, field) || isempty(s.(field))
    return
end
v = s.(field);
if ischar(v) || isstring(v)
    v = cellstr(v);
end
if numel(v) ~= n
    error('saveKcatSource:fields', 'kcatList.%s and kcatList.rxns differ in length.', field)
end
for i = 1:n
    e = v{i};
    if iscell(e)
        e = strjoin(cellfun(@char, e(:)', 'UniformOutput', false), '; ');
    elseif isstring(e)
        e = char(e);
    elseif ~ischar(e)
        e = '';
    end
    col{i} = e;
end
end

function col = numColumn(s, field, n)
col = nan(n, 1);
if isfield(s, field) && ~isempty(s.(field))
    if numel(s.(field)) ~= n
        error('saveKcatSource:fields', 'kcatList.%s and kcatList.rxns differ in length.', field)
    end
    col = double(s.(field)(:));
end
end

function t = formatNumbers(x)
% formatKcatNumber for a vector, in a few vectorized passes.
x = x(:);
t = strings(numel(x), 1);
whole = isfinite(x) & x == round(x) & abs(x) < 1e15;
t(whole) = compose('%d', x(whole));
t(isinf(x) & x > 0) = "Inf";
t(isinf(x) & x < 0) = "-Inf";
open = find(isfinite(x) & ~whole);
for digits = 1:17
    if isempty(open)
        break
    end
    candidate = compose(sprintf('%%.%dg', digits), x(open));
    hit = str2double(candidate) == x(open);
    t(open(hit)) = candidate(hit);
    open = open(~hit);
end
end

function t = cleanText(t)
t = regexprep(char(t), '[\t\r\n]', ' ');
end

function v = geckoVersion()
v = '';
file = fullfile(findGECKOroot, 'version.txt');
if isfile(file)
    v = strtrim(fileread(file));
end
end
