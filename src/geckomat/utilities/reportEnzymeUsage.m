function usageReport = reportEnzymeUsage(ecModel, usageData, highCapUsage)
% reportEnzymeUsage
%   Summarizes the results from enzymeUsage.
%
%  Input:
%   ecModel         a GECKO3 ecModel
%   usageData       output from reportEnzymeUsage
%   highCapUsage    minimum ratio of enzyme capacity usage to be considered
%                   as high usage (Optional, default 0.9)
%
%  Output:
%   usageReport     table with summary information
%
% Usage:
%   usageReport = reportEnzymeUsage(ecModel, usageData, highCapUsage)

if nargin < 3 || isempty(highCapUsage)
    highCapUsage = 0.9
end

usageReport = {};

% Highest usage
highUsageProt = find(usageData.capUsage > highCapUsage);
[~,enzIdx] = ismember(usageData.protID(highUsageProt),ecModel.ec.enzymes);
[row, col] = find(ecModel.ec.rxnEnzMat(:,enzIdx));
[row, ordered] = sort(row);
highUsage.rxnID     = ecModel.ec.rxns(row);
[~, rxnIdx] = ismember(highUsage.rxnID,ecModel.rxns);
highUsage.rxnName   = ecModel.rxnNames(rxnIdx);
protID = highUsageProt(col);
geneID = ecModel.ec.genes(enzIdx(col));
highUsage.protID    = usageData.protID(protID(ordered));
highUsage.geneID    = geneID(ordered);
highUsage.grRules   = ecModel.grRules(rxnIdx);
highUsage.capUsage  = usageData.capUsage(protID(ordered));
highUsage.absUsage  = usageData.absUsage(protID(ordered));

usageReport.highUsage = struct2table(highUsage);

end
