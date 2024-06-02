function mergedKcatList = mergeDLKcatAndFuzzyKcats(kcatListDLKcat, kcatListFuzzy, topOriginLimit, bottomOriginLimit, wildcardLimit)
% mergeDlkcatAndFuzzyKcats
%   Merges the results from DLKcat and fuzzy matching to BRENDA database.
%   Order of preference:
%   1: BRENDA match with correct E.C. number, with origin (see below) not
%      lower than the specified topOriginLimit
%   2: DLKcat match
%   3: BRENDA match with correct E.C. number, with origin below
%      topOriginLimit but not lower than the bottomOriginLimit
%   4: BRENDA match with wildcards in the E.C. number, with not more
%      wildcards than wildcardLimit, and origin not lower than the
%      bottomOriginLimit
%
% Input:
%   kcatListDLKcat      kcatList derived from readDLKcatOutput
%   kcatListFuzzy       kcatList derived from fuzzyKcatMatching
%   topOriginLimit      origin limit for prioritized BRENDA matches. Origin
%                       is explained in more detail below. (Optional,
%                       default 6)
%   bottomOriginLimit   origin limit for low priority BRENDA matches.
%                       Origin is explained in more detail below.
%                       (Optional, default 6)
%   wildcardLimit       maximum number of wildcards in E.C. number of
%                       BRENDA matches (Optional, default 3)
%
% Output:
%   mergedKcatList      merged list of kcats
%   
% The origin parameter:
%   1: correct organism, correct substrate, kcat
%   2: any organism, correct substrate, kcat
%   3: correct organism, any substrate, kcat
%   4: any organism, any substrate, kcat
%   5: correct organism, specific activity
%   6: any organism, specific activity
%
% Example of wildcards in E.C. number:
%   0: 1.1.1.3      glycerol-3-phosphate dehydrogenase (NAD+)
%   1: 1.1.1.-      oxidoreductase, acting on the CH-OH group of donors,
%                   with NAD+ or NADP+ as acceptor
%   2: 1.1.-.-      oxidoreductase, acting on the CH-OH group of donors
%   3: 1.-.-.-      oxidoreductase
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

if (topOriginLimit < 1) || (topOriginLimit > 6)
    error('topPrioOriginLimit should be between 1 and 6.');
end

if (bottomOriginLimit < 1) || (bottomOriginLimit > 6)
    error('originCutLevel should be between 1 and 6.');
end

if (wildcardLimit < 0) || (wildcardLimit > 3)
    error('wildcardCutLevel should be between 0 and 3.');
end

wc = kcatListFuzzy.wildcardLvl;
wc(isnan(wc)) = 1000; %large wildcard

origin = kcatListFuzzy.origin;
origin(isnan(origin)) = 1000; %large origin

prio1 = (wc == 0) & (origin <= topOriginLimit);

rxnsWithPrio1 = unique(kcatListFuzzy.rxns(prio1));

%Things get a bit complicated since not all reactions are in the kcatLists and
%some reactions may appear multiple times
prio2 = true(length(kcatListDLKcat.rxns),1);
prio2(ismember(kcatListDLKcat.rxns, rxnsWithPrio1)) = false;
prio2Rxns = unique(kcatListDLKcat.rxns(prio2));

%The prioritization between wildcards and origin is already done in fuzzy matching,
%so we can join them here
prio3 = ((wc == 0) & (origin > topOriginLimit) & (origin <= bottomOriginLimit)) | ...
        ((wc > 0) & (wc <= wildcardLimit) & (origin <= bottomOriginLimit));
prio3(ismember(kcatListFuzzy.rxns, prio2Rxns)) = false;

fuzzyRxns = prio1 | prio3;

%Now build the merged list, fuzzy followed by dlkcat
%The order of the reactions is therefore not preserved.
mergedKcatList               = struct();
mergedKcatList.source        = 'Merged DLKcat and fuzzy';
[fuzzySrc{1:sum(fuzzyRxns)}] = deal(kcatListFuzzy.source);
[dlkcatSrc{1:sum(prio2)}]    = deal(kcatListDLKcat.source);
mergedKcatList.kcatSource    = [fuzzySrc.';dlkcatSrc.'];
mergedKcatList.rxns          = [kcatListFuzzy.rxns(fuzzyRxns);kcatListDLKcat.rxns(prio2)];
mergedKcatList.genes         = [cell(sum(fuzzyRxns),1);kcatListDLKcat.genes(prio2)];
mergedKcatList.substrates    = [kcatListFuzzy.substrates(fuzzyRxns);kcatListDLKcat.substrates(prio2)];
mergedKcatList.kcats         = [kcatListFuzzy.kcats(fuzzyRxns);kcatListDLKcat.kcats(prio2)];
mergedKcatList.eccodes       = [kcatListFuzzy.eccodes(fuzzyRxns);cell(sum(prio2),1)];
mergedKcatList.wildcardLvl   = [kcatListFuzzy.wildcardLvl(fuzzyRxns);nan(sum(prio2),1)];
mergedKcatList.origin        = [kcatListFuzzy.origin(fuzzyRxns);nan(sum(prio2),1)];
end
