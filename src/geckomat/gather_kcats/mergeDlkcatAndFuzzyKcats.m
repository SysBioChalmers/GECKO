function mergedKcatList = mergeDlkcatAndFuzzyKcats(kcatListDlkcat, kcatListFuzzy, topPrioOriginLimit, originCutLevel, wildcardCutLevel)
% mergeDlkcatAndFuzzyKcats
%   Merges the results from dlKcat and fuzzy matching, where the direct matches
%   from fuzzy matching down to an origin value of topPrioOriginLimit has the 
%   highest priority, followed by dlKcat, and then by the fuzzy direct matches below
%   topPrioOriginLimit, but above or equal to originCutLevel. After that fuzzy 
%   wildcard matches with wildcard >= wildcardCutLevel follows.
%
% Input:
%   kcatListDlkcat              kcats from dlkcat
%   kcatListFuzzy               kcats from fuzzy matching
%   topPrioOriginLimit          described above, optional, default 6
%   originCutLevel              described above, optional, default 6
%   wildcardCutLevel            described above, default 0 (discards all wildcards)
%
% Output:
%   mergedKcatList  The merged list of kcats
%   

% The origin parameter:
%   1: correct organism, correct substrate, kcat
%   2: any organism, correct substrate, kcat
%   3: correct organism, any substrate, kcat
%   4: any organism, any substrate, kcat
%   5: correct organism, specific activity
%   6: any organism, specific activity


if nargin < 5
    wildcardCutLevel = 0;
end

if nargin < 4
    originCutLevel = 6;
end

if nargin < 3
    topPrioOriginLimit = 6;
end

if (topPrioOriginLimit < 1) || (topPrioOriginLimit > 6)
    error('topPrioOriginLimit should be between 1 and 6.');
end

if (originCutLevel < 1) || (originCutLevel > 6)
    error('originCutLevel should be between 1 and 6.');
end

if (wildcardCutLevel < 0) || (wildcardCutLevel > 3)
    error('wildcardCutLevel should be between 0 and 3.');
end

wc = kcatListFuzzy.wildcardLvl;
wc(isnan(wc)) = 1000; %large wildcard

origin = kcatListFuzzy.origin;
origin(isnan(origin)) = 1000; %large origin

prio1 = (wc == 0) & (origin <= topPrioOriginLimit);

rxnsWithPrio1 = unique(kcatListFuzzy.rxns(prio1));

%Things get a bit complicated since not all reactions are in the kcatLists and
%some reactions may appear multiple times
prio2 = true(length(kcatListDlkcat.rxns),1);
prio2(ismember(kcatListDlkcat.rxns, rxnsWithPrio1)) = false;
prio2Rxns = unique(kcatListDlkcat.rxns(prio2));

%The prioritization between wildcards and origin is already done in fuzzy matching,
%so we can join them here
prio3 = ((wc == 0) & (origin > topPrioOriginLimit) & (origin <= originCutLevel)) | ...
        ((wc > 0) & (wc <= wildcardCutLevel) & (origin <= originCutLevel));
prio3(ismember(kcatListFuzzy.rxns, prio2Rxns)) = false;

fuzzyRxns = prio1 | prio3;

%Now build the merged list, fuzzy followed by dlkcat
%The order of the reactions is therefore not preserved.
mergedKcatList = struct();
mergedKcatList.source = 'Merged dlkcat and fuzzy';
[fuzzySrc{1:sum(fuzzyRxns)}] = deal(kcatListFuzzy.source);
[dlkcatSrc{1:sum(prio2)}] = deal(kcatListDlkcat.source);
mergedKcatList.kcatSource = [fuzzySrc.';dlkcatSrc.'];
mergedKcatList.rxns = [kcatListFuzzy.rxns(fuzzyRxns);kcatListDlkcat.rxns(prio2)];
mergedKcatList.genes = [cell(sum(fuzzyRxns),1);kcatListDlkcat.genes(prio2)];
mergedKcatList.substrates = [kcatListFuzzy.substrates(fuzzyRxns);kcatListDlkcat.substrates(prio2)];
mergedKcatList.kcats = [kcatListFuzzy.kcats(fuzzyRxns);kcatListDlkcat.kcats(prio2)];
mergedKcatList.eccodes = [kcatListFuzzy.eccodes(fuzzyRxns);cell(sum(prio2),1)];
mergedKcatList.wildcardLvl = [kcatListFuzzy.wildcardLvl(fuzzyRxns);nan(sum(prio2),1)];
mergedKcatList.origin = [kcatListFuzzy.origin(fuzzyRxns);nan(sum(prio2),1)];

end
