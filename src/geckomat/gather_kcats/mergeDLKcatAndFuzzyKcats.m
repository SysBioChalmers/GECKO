function mergedKcatList = mergeDLKcatAndFuzzyKcats(kcatListDLKcat, kcatListFuzzy, topOriginLimit, bottomOriginLimit, wildcardLimit)
% mergeDLKcatAndFuzzyKcats
%   Deprecated thin wrapper around mergeKcats, kept for the common
%   DLKcat + fuzzy-BRENDA case. New code should call mergeKcats directly,
%   which also handles a single list that mixes several sources (e.g. an
%   OpenKineticsPredictor result with BRENDA/Sabio-RK/CataPro rows).
%
%   Merges the results from DLKcat and fuzzy matching to BRENDA database.
%   Order of preference (best first):
%   1: BRENDA match with correct E.C. number and origin not lower than
%      topOriginLimit            (tier 'database_top')
%   2: DLKcat match              (tier 'dlkcat')
%   3: BRENDA match with correct E.C. number, origin below topOriginLimit
%      but not lower than bottomOriginLimit, or with no more wildcards
%      than wildcardLimit        (tier 'database_bottom')
%
% Input:
%   kcatListDLKcat      kcatList derived from readDLKcatOutput
%   kcatListFuzzy       kcatList derived from fuzzyKcatMatching
%   topOriginLimit      origin limit for prioritized BRENDA matches
%                       (Optional, default 6)
%   bottomOriginLimit   origin limit for low priority BRENDA matches
%                       (Optional, default 6)
%   wildcardLimit       maximum number of wildcards in E.C. number of
%                       BRENDA matches (Optional, default 3)
%
% Output:
%   mergedKcatList      merged list of kcats
%
% Usage:
%   mergedKcatList = mergeDLKcatAndFuzzyKcats(kcatListDLKcat, kcatListFuzzy, topOriginLimit, bottomOriginLimit, wildcardLimit)

if nargin < 5
    wildcardLimit = 3;
end
if nargin < 4
    bottomOriginLimit = 6;
end
if nargin < 3
    topOriginLimit = 6;
end

warning('mergeDLKcatAndFuzzyKcats:deprecated', ...
    ['mergeDLKcatAndFuzzyKcats is deprecated; use mergeKcats instead. ' ...
     'The old name will be removed in a future release.']);

% Fuzzy list first, then DLKcat, so surviving rows keep the legacy order
% (fuzzy block followed by DLKcat block).
mergedKcatList = mergeKcats({kcatListFuzzy, kcatListDLKcat}, ...
    {'database_top', 'dlkcat', 'database_bottom'}, ...
    topOriginLimit, bottomOriginLimit, wildcardLimit);
end
