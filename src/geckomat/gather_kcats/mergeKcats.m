function mergedKcatList = mergeKcats(kcatLists, sourcePriority, topOriginLimit, bottomOriginLimit, wildcardLimit, databaseSources)
% mergeKcats
%   Merges any number of kcatLists, keeping the highest-priority source
%   for each reaction. Generalizes mergeDLKcatAndFuzzyKcats: instead of
%   two fixed inputs (DLKcat + fuzzy BRENDA) it accepts a cell array of
%   kcatLists, each of which may itself mix several sources (for example
%   the output of fetchOpenKineticsPredictor, carrying BRENDA, Sabio-RK
%   and CataPro values in its kcatSource field). A caller-supplied
%   sourcePriority decides which source wins for each reaction.
%
% Input:
%   kcatLists           cell array of kcatList structs (as produced by
%                       fuzzyKcatMatching, readDLKcatOutput,
%                       fetchOpenKineticsPredictor, ...)
%   sourcePriority      cell array of tier tokens / source labels, best
%                       first. Three reserved tier tokens may be used
%                       instead of a literal source label:
%                         'database_exact'  an exact experimental-database
%                                           measurement (BRENDA/Sabio-RK
%                                           row with no fuzzy metadata, i.e.
%                                           a direct OKP hit, no wildcards).
%                         'database_top'    a fuzzy BRENDA match with
%                                           wildcardLvl==0 and
%                                           origin<=topOriginLimit.
%                         'database_bottom' weaker fuzzy BRENDA matches
%                                           within wildcardLimit /
%                                           bottomOriginLimit.
%                       Any other token is matched against each row's
%                       source after folding to lowercase (e.g. 'dlkcat',
%                       'catapro'). Sources not listed are dropped.
%   topOriginLimit      origin limit for database_top (Optional, default 6)
%   bottomOriginLimit   origin limit for database_bottom (Optional, def 6)
%   wildcardLimit       maximum wildcards for database_bottom's wildcard
%                       branch (Optional, default 3)
%   databaseSources     source labels treated as experimental databases,
%                       i.e. routed to the database_* tiers (Optional,
%                       default {'brenda','sabio_rk'})
%
% Output:
%   mergedKcatList      merged kcatList. Per-row provenance is preserved
%                       verbatim in the kcatSource field; the scalar
%                       source field is set to 'merged'. Row order follows
%                       the order of the input lists (each list's surviving
%                       rows, in turn), so a two-list call behaves like the
%                       legacy mergeDLKcatAndFuzzyKcats.
%
% The origin parameter:
%   1: correct organism, correct substrate, kcat
%   2: any organism, correct substrate, kcat
%   3: correct organism, any substrate, kcat
%   4: any organism, any substrate, kcat
%   5: correct organism, specific activity
%   6: any organism, specific activity
%
% Usage:
%   mergedKcatList = mergeKcats(kcatLists, sourcePriority, topOriginLimit, bottomOriginLimit, wildcardLimit, databaseSources)

if nargin < 6 || isempty(databaseSources)
    databaseSources = {'brenda','sabio_rk'};
end
if nargin < 5 || isempty(wildcardLimit)
    wildcardLimit = 3;
end
if nargin < 4 || isempty(bottomOriginLimit)
    bottomOriginLimit = 6;
end
if nargin < 3 || isempty(topOriginLimit)
    topOriginLimit = 6;
end

if (topOriginLimit < 1) || (topOriginLimit > 6)
    error('topOriginLimit should be between 1 and 6.');
end
if (bottomOriginLimit < 1) || (bottomOriginLimit > 6)
    error('bottomOriginLimit should be between 1 and 6.');
end
if (wildcardLimit < 0) || (wildcardLimit > 3)
    error('wildcardLimit should be between 0 and 3.');
end
if isempty(sourcePriority)
    error('sourcePriority must be a non-empty cell array.');
end
if ~iscell(kcatLists)
    error('kcatLists must be a cell array of kcatList structs.');
end

dbSources = cellfun(@normalizeSource, databaseSources, 'UniformOutput', false);

% Rank map: normalised priority token -> rank (1 = best). First wins.
rankMap = containers.Map('KeyType','char','ValueType','double');
for r = 1:numel(sourcePriority)
    key = normalizeSource(sourcePriority{r});
    if ~isKey(rankMap, key)
        rankMap(key) = r;
    end
end

% Concatenate all lists, normalising missing fields per row.
rxns = {}; kcats = []; genes = {}; substrates = {};
eccodes = {}; wcl = []; origin = []; rawSource = {};
for k = 1:numel(kcatLists)
    s = kcatLists{k};
    if isempty(s) || ~isfield(s,'rxns') || isempty(s.rxns)
        continue
    end
    n          = numel(s.rxns);
    rxns       = [rxns;       s.rxns(:)];
    kcats      = [kcats;      double(s.kcats(:))];
    genes      = [genes;      getColCell(s, 'genes', n)];
    substrates = [substrates; getColCell(s, 'substrates', n)];
    eccodes    = [eccodes;    getColCell(s, 'eccodes', n)];
    wcl        = [wcl;        getColNum(s, 'wildcardLvl', n)];
    origin     = [origin;     getColNum(s, 'origin', n)];
    rawSource  = [rawSource;  getRowSource(s, n)];
end

N = numel(rxns);
if N == 0
    mergedKcatList = emptyMergedList();
    return
end

% Drop rows without a usable turnover number. This also removes the
% kcats==0 rows fuzzyKcatMatching emits for unmatched reactions, so a
% remaining NaN-metadata database row is unambiguously an exact OKP hit.
valid          = kcats > 0;
rxns           = rxns(valid);
kcats          = kcats(valid);
genes          = genes(valid);
substrates     = substrates(valid);
eccodes        = eccodes(valid);
wcl            = wcl(valid);
origin         = origin(valid);
rawSource      = rawSource(valid);
N              = numel(rxns);
if N == 0
    mergedKcatList = emptyMergedList();
    return
end

token  = cellfun(@normalizeSource, rawSource, 'UniformOutput', false);
isDb   = ismember(token, dbSources);
hasMeta = ~isnan(wcl) & ~isnan(origin);

% Tier per row ('' marks a database row that fails the quality gate).
tier = token;
dbIdx = find(isDb);
for ii = dbIdx'
    if ~hasMeta(ii)
        tier{ii} = 'database_exact';
    elseif wcl(ii) == 0 && origin(ii) <= topOriginLimit
        tier{ii} = 'database_top';
    elseif (wcl(ii) == 0 && origin(ii) > topOriginLimit && origin(ii) <= bottomOriginLimit) || ...
           (wcl(ii) > 0  && wcl(ii) <= wildcardLimit       && origin(ii) <= bottomOriginLimit)
        tier{ii} = 'database_bottom';
    else
        tier{ii} = '';
    end
end

% Rank per row (NaN if the tier is not listed in sourcePriority).
rankPerRow = nan(N,1);
for ii = 1:N
    if ~isempty(tier{ii}) && isKey(rankMap, tier{ii})
        rankPerRow(ii) = rankMap(tier{ii});
    end
end

% Warn about non-database sources that were dropped for not being listed.
droppedTokens = unique(token(~isDb & isnan(rankPerRow)));
if ~isempty(droppedTokens)
    warning('mergeKcats:droppedSources', ...
        'Dropping rows from source(s) not listed in sourcePriority: %s', ...
        strjoin(droppedTokens', ', '));
end

% For each reaction keep all rows of its best (lowest-rank) tier.
keep = false(N,1);
validRank = ~isnan(rankPerRow);
uRxns = unique(rxns(validRank));
for k = 1:numel(uRxns)
    ri = find(strcmp(rxns, uRxns{k}) & validRank);
    minR = min(rankPerRow(ri));
    keep(ri(rankPerRow(ri) == minR)) = true;
end

mergedKcatList             = struct();
mergedKcatList.source      = 'merged';
mergedKcatList.kcatSource  = rawSource(keep);
mergedKcatList.rxns        = rxns(keep);
mergedKcatList.genes       = genes(keep);
mergedKcatList.substrates  = substrates(keep);
mergedKcatList.kcats       = kcats(keep);
mergedKcatList.eccodes     = eccodes(keep);
mergedKcatList.wildcardLvl = wcl(keep);
mergedKcatList.origin      = origin(keep);
end

% ------------------------------------------------------------------------- %
function s = normalizeSource(label)
% Fold a source label to lowercase snake_case ('Sabio-RK' -> 'sabio_rk').
s = lower(strtrim(char(label)));
s = regexprep(s, '[^0-9a-z]+', '_');
s = regexprep(s, '^_+|_+$', '');
end

function v = getColCell(s, name, n)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name)(:);
    if numel(v) == 1 && n > 1
        v = repmat(v, n, 1);
    end
else
    v = repmat({[]}, n, 1);
end
end

function v = getColNum(s, name, n)
if isfield(s, name) && ~isempty(s.(name))
    v = double(s.(name)(:));
else
    v = nan(n, 1);
end
end

function v = getRowSource(s, n)
if isfield(s, 'kcatSource') && ~isempty(s.kcatSource)
    v = s.kcatSource(:);
elseif isfield(s, 'source')
    v = repmat({char(s.source)}, n, 1);
else
    v = repmat({''}, n, 1);
end
end

function mergedKcatList = emptyMergedList()
mergedKcatList             = struct();
mergedKcatList.source      = 'merged';
mergedKcatList.kcatSource  = {};
mergedKcatList.rxns        = {};
mergedKcatList.genes       = {};
mergedKcatList.substrates  = {};
mergedKcatList.kcats       = [];
mergedKcatList.eccodes     = {};
mergedKcatList.wildcardLvl = [];
mergedKcatList.origin      = [];
end
