function findings = reviewKcatAssignment(rxnIds,kcats,sources,ecCodes,varargin)
% reviewKcatAssignment  Flag kcat assignments that deserve a second look.
%
% A tuner's largest corrections are rarely discoveries about biology.
% They are places the assignment went wrong: a value derived from one
% specific-activity measurement, a database maximum standing in for a
% distribution with a two-hundred-fold spread, one fallback number reused
% across a whole family of reactions. Finding those directly is cheaper
% than inferring them from a fit, and a corrected prior propagates to
% every model built from the same databases, where a tuned kcat stays in
% the model it was tuned on.
%
% Two design choices, both from measurement rather than taste (see
% geckopy's gather_kcats/review.py, which this ports):
%
% Rank by leverage, and cut the report there. Tightening the checks
% themselves does not shorten a report nearly as well, because looking
% odd and mattering are unrelated properties. Filtering by leverage keeps
% the list short and keeps the entries worth a curator's time.
%
% Report a parameter once. Isozyme copies the assignment could not
% separate share a value and a source, so they share a finding; listing
% them separately fills the report with duplicates of one decision.
%
% Parameters
% ----------
% rxnIds, kcats, sources, ecCodes : cell/double
%     ec.rxns, ec.kcat, ec.source, ec.eccodes (or the equivalent
%     subset), all the same length.
%
% Name-Value Arguments
% --------------------
% ecStats : containers.Map
%     EC code -> struct(nKcat,kcatMax,kcatMedian,values,nSa), as returned
%     by ecStatsFromBrenda. Checks needing it are skipped when it is
%     absent, so this still runs without a database.
% enzymes : cell
%     protein ids per kcat (cell of cell arrays), used to count how many
%     other reactions share an enzyme.
% names : cell
%     reaction names, same length as rxnIds (default: not shown).
% leverage : double
%     per-kcat effect on the distance from a one-at-a-time screen (e.g.
%     screenKcatLeverage's leverage column, expanded back to
%     ecModel.ec.rxns positions). Without it findings cannot be ranked
%     and are returned in input order.
% slow, fast : double
%     magnitude check bounds (default 1e-2, 1e4).
% ecMaxGap : double
%     how many times the EC median a value must exceed, while sitting
%     exactly at the EC max, to flag as ec-maximum (default 3.0).
% repeatMin : double
%     how many distinct reactions must share one kcat value to flag as
%     repeated-value (default 5).
% include : cell
%     names from {'standard-fallback','tied-isozymes'} to report as well
%     (default: neither -- both are deliberate choices, not mistakes, so
%     stay quiet unless asked for).
% coverage : double
%     truncate the report to the rows carrying this share of the
%     *flagged* leverage (default: no truncation).
% top : double
%     truncate the report to this many rows (default: no truncation).
%
% Returns
% -------
% findings : table
%     one row per free parameter: rxnId, name, ecCode, source, kcat,
%     leverage, nIsozymes, otherReactions, checks (comma-joined string),
%     detail. Write to TSV with e.g.
%         writetable(findings,'findings.tsv','FileType','text','Delimiter','\t')
%
% See also
% --------
% ecStatsFromBrenda, isozymeTieMap, screenKcatLeverage
p = parseGECKOargs(varargin, { ...
    'ecStats', containers.Map('KeyType','char','ValueType','any'); ...
    'enzymes', {}; ...
    'names', {}; ...
    'leverage', []; ...
    'slow', 1e-2; ...
    'fast', 1e4; ...
    'ecMaxGap', 3.0; ...
    'repeatMin', 5; ...
    'include', {}; ...
    'coverage', []; ...
    'top', []});
ecStats   = p.ecStats;
enzymes   = p.enzymes;
names     = p.names;
leverage  = p.leverage;
slow      = p.slow;
fast      = p.fast;
ecMaxGap  = p.ecMaxGap;
repeatMin = p.repeatMin;
include   = p.include;
coverage  = p.coverage;
top       = p.top;

kcats = kcats(:);
n = numel(kcats);
if numel(rxnIds) ~= n || numel(sources) ~= n || numel(ecCodes) ~= n
    error('rxnIds, sources and ecCodes must all have %d entries, matching kcats.', n)
end
if isempty(leverage)
    lev = zeros(n,1);
else
    lev = leverage(:);
end

tie = isozymeTieMap(rxnIds,kcats,'sources',sources);
reps = unique(tie);
memberLists = cell(numel(reps),1);
for r = 1:numel(reps)
    memberLists{r} = find(tie == reps(r));
end

% A value reused across many distinct reactions is a fallback wearing a
% measurement's label, whatever its source says.
reuseMap = containers.Map('KeyType','char','ValueType','any');
for i = 1:n
    key = reuseKey(kcats(i));
    base = baseReactionLocal(rxnIds{i});
    if isKey(reuseMap,key)
        reuseMap(key) = union(reuseMap(key), {base});
    else
        reuseMap(key) = {base};
    end
end

enzymeRxns = containers.Map('KeyType','char','ValueType','any');
if ~isempty(enzymes)
    for i = 1:n
        for e = 1:numel(enzymes{i})
            prot = char(enzymes{i}{e});
            base = baseReactionLocal(rxnIds{i});
            if isKey(enzymeRxns,prot)
                enzymeRxns(prot) = union(enzymeRxns(prot), {base});
            else
                enzymeRxns(prot) = {base};
            end
        end
    end
end

rxnIdOut=[]; nameOut=[]; ecOut=[]; sourceOut=[]; kcatOut=[]; leverageOut=[];
nIsoOut=[]; otherOut=[]; checksOut=[]; detailOut=[];

for r = 1:numel(reps)
    idx = memberLists{r};
    i = idx(1);
    ec = firstEc(ecCodes{i});
    st = [];
    if isKey(ecStats,ec)
        st = ecStats(ec);
    end
    checks = {}; detail = {};

    if ~isempty(st) && st.nKcat == 0
        checks{end+1} = 'no-ec-evidence'; %#ok<AGROW>
        if st.nSa > 0
            detail{end+1} = sprintf('no kcat rows for EC %s; %d specific-activity rows instead', ec, st.nSa); %#ok<AGROW>
        else
            detail{end+1} = sprintf('no kcat rows for EC %s', ec); %#ok<AGROW>
        end
    end
    key = reuseKey(kcats(i));
    if isKey(reuseMap,key) && numel(reuseMap(key)) >= repeatMin
        checks{end+1} = 'repeated-value'; %#ok<AGROW>
        detail{end+1} = sprintf('%g 1/s reused across %d reactions', kcats(i), numel(reuseMap(key))); %#ok<AGROW>
    end
    if ~isempty(st) && st.nKcat > 0 && isfinite(st.kcatMedian) && st.kcatMedian > 0 ...
            && abs(kcats(i)-st.kcatMax) <= 1e-9*max(kcats(i),1) && st.kcatMax/st.kcatMedian > ecMaxGap
        checks{end+1} = 'ec-maximum'; %#ok<AGROW>
        detail{end+1} = sprintf('took the EC maximum, %.0fx its median over %d rows', st.kcatMax/st.kcatMedian, st.nKcat); %#ok<AGROW>
    end
    if kcats(i) < slow || kcats(i) > fast
        checks{end+1} = 'magnitude'; %#ok<AGROW>
        detail{end+1} = sprintf('%g 1/s is outside %g-%g', kcats(i), slow, fast); %#ok<AGROW>
    end
    if strcmpi(sources{i},'custom') && ~isempty(st) && ismember(round(kcats(i),6), st.values)
        checks{end+1} = 'custom-duplicate'; %#ok<AGROW>
        detail{end+1} = 'custom value equals a database value for its own EC'; %#ok<AGROW>
    end
    if numel(idx) > 1
        if ismember('tied-isozymes',include)
            checks{end+1} = 'tied-isozymes'; %#ok<AGROW>
        end
        detail{end+1} = sprintf('%d isozyme copies share this value', numel(idx)); %#ok<AGROW>
    end
    if strcmpi(sources{i},'standard')
        if ~ismember('standard-fallback',include)
            continue
        end
        checks{end+1} = 'standard-fallback'; %#ok<AGROW>
        detail{end+1} = 'the model''s fallback value, not a measurement'; %#ok<AGROW>
    end

    if isempty(checks)
        continue
    end
    others = {};
    if ~isempty(enzymes)
        for e = 1:numel(enzymes{i})
            prot = char(enzymes{i}{e});
            if isKey(enzymeRxns,prot)
                others = union(others, enzymeRxns(prot));
            end
        end
    end
    others = setdiff(others, {baseReactionLocal(rxnIds{i})});

    rxnIdOut{end+1,1}   = rxnIds{i}; %#ok<AGROW>
    if isempty(names)
        nameOut{end+1,1} = ''; %#ok<AGROW>
    else
        nameOut{end+1,1} = names{i}; %#ok<AGROW>
    end
    ecOut{end+1,1}       = ec; %#ok<AGROW>
    sourceOut{end+1,1}   = sources{i}; %#ok<AGROW>
    kcatOut(end+1,1)     = kcats(i); %#ok<AGROW>
    leverageOut(end+1,1) = max(lev(idx)); %#ok<AGROW>
    nIsoOut(end+1,1)     = numel(idx); %#ok<AGROW>
    otherOut(end+1,1)    = numel(others); %#ok<AGROW>
    checksOut{end+1,1}   = strjoin(checks,','); %#ok<AGROW>
    detailOut{end+1,1}   = strjoin(detail,'; '); %#ok<AGROW>
end

if isempty(rxnIdOut)
    findings = table(cell(0,1),cell(0,1),cell(0,1),cell(0,1),zeros(0,1),zeros(0,1), ...
        zeros(0,1),zeros(0,1),cell(0,1),cell(0,1), 'VariableNames', ...
        {'rxnId','name','ecCode','source','kcat','leverage','nIsozymes','otherReactions','checks','detail'});
    return
end

findings = table(rxnIdOut,nameOut,ecOut,sourceOut,kcatOut,leverageOut,nIsoOut,otherOut,checksOut,detailOut, ...
    'VariableNames',{'rxnId','name','ecCode','source','kcat','leverage','nIsozymes','otherReactions','checks','detail'});
findings = sortrows(findings,'leverage','descend');

flaggedTotal = sum(findings.leverage);
if ~isempty(coverage) && flaggedTotal > 0
    acc = cumsum(findings.leverage);
    keepN = find(acc/flaggedTotal >= coverage, 1, 'first');
    if isempty(keepN)
        keepN = height(findings);
    end
    findings = findings(1:keepN,:);
end
if ~isempty(top)
    findings = findings(1:min(top,height(findings)),:);
end
end

function key = reuseKey(kcat)
if kcat > 0
    key = sprintf('%.9f', round(log(kcat),9));
else
    key = 'NONE';
end
end

function ec = firstEc(ecCode)
parts = strsplit(char(ecCode), ';');
ec = strtrim(parts{1});
end

function base = baseReactionLocal(rxnId)
% Same one-line split as isozymeTieMap's local baseReaction. Kept as a
% local duplicate rather than exported from isozymeTieMap.m, since
% MATLAB only exposes a file's first function to other files and this is
% its only caller here.
parts = strsplit(rxnId, '_EXP_');
base = parts{1};
end
