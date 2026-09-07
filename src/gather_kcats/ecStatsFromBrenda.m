function ecStats = ecStatsFromBrenda(KCATcell,SAcell)
% ecStatsFromBrenda  Summarise loadBRENDAdata's output per EC number, for
% reviewKcatAssignment's ecStats argument.
%
% TODO: verify against a real BRENDA pull before relying on this. Ported
% from geckopy's ec_stats_from_brenda without being able to run it
% against real BRENDA data in this environment. In particular,
% loadBRENDAdata's KCATcell only carries the BRENDA *max* aggregate per
% (EC, substrate, organism) triple (column 4 of kcat.tsv) -- the median
% column is read internally but discarded before KCATcell is returned --
% so kcatMedian below is always NaN, and reviewKcatAssignment's
% "ec-maximum" check (which needs a finite kcatMedian) will never fire
% until loadBRENDAdata is extended to also expose it.
%
% Parameters
% ----------
% KCATcell : cell
%     as returned by loadBRENDAdata's first output: {ecCode, substrate,
%     organism, kcatMax}.
% SAcell : cell
%     as returned by loadBRENDAdata's second output: {ecCode, organism,
%     kcatFromSA, mw}.
%
% Returns
% -------
% ecStats : containers.Map
%     EC code -> struct(nKcat,kcatMax,kcatMedian,values,nSa), keyed by
%     the same (uppercased) EC strings BRENDA uses.
%
% See also
% --------
% reviewKcatAssignment, loadBRENDAdata
ecCol  = upper(KCATcell{1});
maxCol = KCATcell{4};

saEcCol = {};
if ~isempty(SAcell) && numel(SAcell) >= 1 && ~isempty(SAcell{1})
    saEcCol = upper(SAcell{1});
end

ecList = unique([ecCol(:); saEcCol(:)]);
ecStats = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(ecList)
    ec = char(ecList{i});
    rows = strcmp(ecCol,ec);
    nKcat = sum(rows);
    if nKcat > 0
        kcatMax = max(maxCol(rows));
        values  = unique(round(maxCol(rows),6));
    else
        kcatMax = NaN;
        values  = [];
    end
    nSa = sum(strcmp(saEcCol,ec));
    ecStats(ec) = struct('nKcat',nKcat,'kcatMax',kcatMax,'kcatMedian',NaN,'values',values,'nSa',nSa);
end
end
