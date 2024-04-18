function usageReport = topEnzymeUsage(ecModel, usageData, highCapUsage, topAbsUsage)
% reportEnzymeUsage
%   Summarizes the results from enzymeUsage.
%
%  Input:
%   ecModel         a GECKO3 ecModel
%   usageData       output from reportEnzymeUsage
%   highCapUsage    minimum ratio of enzyme capacity usage to be considered
%                   as high usage (Optional, default 0.9, refering to a 
%                   minimum of 90% capacity usage)
%   topAbsUsage     number of top enzymes with highest absolute usage
%                   (Optional, default 10, returning the top 10 enzymes
%                   with highest absolute usage. With Inf or 0, all enzymes
%                   are returned)
%
%  Output:
%   usageReport     table with summary information
%
% Usage:
%   usageReport = reportEnzymeUsage(ecModel, usageData, highCapUsage)

if nargin < 3 || isempty(highCapUsage)
    highCapUsage = 0.9;
end
if nargin < 4 || isempty(topAbsUsage)
    topAbsUsage = 10;
end

usageReport = {};

% Highest capacity usage
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

usageReport.highCapUsage = struct2table(highUsage);

% Highest absolute usage
[~,topUse]          = sort(usageData.absUsage,'descend');
topEnzyme           = usageData.protID(topUse(1:topAbsUsage));
[~,b]       = ismember(topEnzyme,ecModel.ec.enzymes);
geneIDs     = ecModel.ec.genes(b);
topUsage.protID     = {};
topUsage.geneID     = {};
topUsage.absUsage   = [];
topUsage.percUsage  = [];
topUsage.kcat       = [];
topUsage.source     = {};
topUsage.rxnID      = {};
topUsage.rxnNames   = {};
topUsage.grRules    = {};

% Calculate the protein pool flux from the 'prot_pool_exchange' reaction
protPoolExchangeFlux = -ecModel.lb(strcmp(ecModel.rxns,'prot_pool_exchange'));

% Obtain the fluxes
sol = solveLP(ecModel);
if isempty(sol.x)
    error('No solution found.');
end
fluxValues = sol.x;

% Sum fluxes for all 'usage_prot_' reactions, excluding the 'usage_prot_standard'
usageProtIndices = startsWith(ecModel.rxns, 'usage_prot_') & ...
                   ~contains(ecModel.rxns, 'standard');

% Sum the absolute values of the usage fluxes
totalUsageProtFlux = sum(abs(fluxValues(usageProtIndices)));

% Define the new protein pool as the sum of prot_pool_exchange flux and total usage_prot fluxes
protPool = (protPoolExchangeFlux + totalUsageProtFlux)/100;

fprintf('Total Protein Pool Flux: %f\n', protPoolExchangeFlux);
fprintf('total UsageProt Flux: %f\n', totalUsageProtFlux);


for i=1:numel(topEnzyme)
    [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel,topEnzyme{i});
    rxnNumber = numel(rxns);
    topUsage.protID(end+1:end+rxnNumber,1)      = topEnzyme(i);
    topUsage.geneID(end+1:end+rxnNumber,1)      = geneIDs(i);
    topUsage.absUsage(end+1:end+rxnNumber,1)    = usageData.absUsage(topUse(i));
    topUsage.percUsage(end+1:end+rxnNumber,1)   = topUsage.absUsage(end-(rxnNumber-1):end,1)/protPool;
    topUsage.kcat(end+1:end+rxnNumber,1)        = kcat;
    topUsage.source(end+1:end+rxnNumber,1)      = ecModel.ec.source(idx);
    topUsage.rxnID(end+1:end+rxnNumber,1)       = rxns;
    topUsage.rxnNames(end+1:end+rxnNumber,1)    = rxnNames;
    topUsage.grRules(end+1:end+rxnNumber,1)     = grRules;
end
usageReport.topAbsUsage = struct2table(topUsage);
end
