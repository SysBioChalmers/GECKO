function usageReport = reportEnzymeUsage(ecModel, usageData, highCapUsage, topAbsUsage)
% reportEnzymeUsage
%   Summarizes the results from enzymeUsage.
%
%  Input:
%   ecModel         a GECKO3 ecModel
%   usageData       output from enzymeUsage
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
highEnzyme    = usageData.protID(highUsageProt);
[~,enzIdx]    = ismember(highEnzyme,ecModel.ec.enzymes);
geneIDs       = ecModel.ec.genes(enzIdx);

highUsage.protID     = {};
highUsage.geneID     = {};
highUsage.absUsage   = [];
highUsage.capUsage   = [];
highUsage.kcat       = [];
highUsage.source     = {};
highUsage.rxnID      = {};
highUsage.rxnNames   = {};
highUsage.grRules    = {};

for i=1:numel(enzIdx)
    [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel,ecModel.ec.enzymes(enzIdx(i)));
    % See if all reactions carried flux
    [~,rIdx] = ismember(rxns,ecModel.rxns);
    carriedFlux = usageData.fluxes(rIdx) > 1e-7;
    if isscalar(find(carriedFlux))
        highUsage.protID(end+1,1)      = highEnzyme(i);
        highUsage.geneID(end+1,1)      = geneIDs(i);
        highUsage.absUsage(end+1,1)    = usageData.absUsage(enzIdx(i));
        highUsage.capUsage(end+1,1)    = usageData.capUsage(enzIdx(i));
        highUsage.kcat(end+1,1)        = kcat(carriedFlux);
        highUsage.source(end+1,1)      = ecModel.ec.source(idx(carriedFlux));
        highUsage.rxnID(end+1,1)       = rxns(carriedFlux);
        highUsage.rxnNames(end+1,1)    = rxnNames(carriedFlux);
        highUsage.grRules(end+1,1)     = grRules(carriedFlux);
    else
        % Add one entry for combined usage
        highUsage.protID(end+1,1)      = highEnzyme(i);
        highUsage.geneID(end+1,1)      = geneIDs(i);
        highUsage.absUsage(end+1,1)    = usageData.absUsage(enzIdx(i));
        highUsage.capUsage(end+1,1)    = usageData.capUsage(enzIdx(i));
        highUsage.kcat(end+1,1)        = nan;
        highUsage.source{end+1,1}      = '===';
        highUsage.rxnID{end+1,1}       = '===';
        highUsage.rxnNames{end+1,1}    = 'involved in multiple rxns, usage combined, individual rxns below';
        highUsage.grRules{end+1,1}     = '===';
        % Recalculate reaction-specific usage

        rIdx = rIdx(carriedFlux);
        enzFlux = usageData.fluxes(rIdx);
        enzMet = strcat('prot_',highEnzyme{i});
        [~, enzEcIdx] = ismember(enzMet,ecModel.mets);
        indAbsUse = full(transpose(-ecModel.S(enzEcIdx,rIdx)).*enzFlux);
        indCapUse = (indAbsUse /sum(indAbsUse)) * usageData.capUsage(enzIdx(i));

        rxnNumber = length(rIdx);
        highUsage.protID(end+1:end+rxnNumber,1)      = highEnzyme(i);
        highUsage.geneID(end+1:end+rxnNumber,1)      = geneIDs(i);
        highUsage.absUsage(end+1:end+rxnNumber,1)    = indAbsUse;
        highUsage.capUsage(end+1:end+rxnNumber,1)    = indCapUse;
        highUsage.kcat(end+1:end+rxnNumber,1)        = kcat(carriedFlux);
        highUsage.source(end+1:end+rxnNumber,1)      = ecModel.ec.source(idx(carriedFlux));
        highUsage.rxnID(end+1:end+rxnNumber,1)       = rxns(carriedFlux);
        highUsage.rxnNames(end+1:end+rxnNumber,1)    = rxnNames(carriedFlux);
        highUsage.grRules(end+1:end+rxnNumber,1)     = grRules(carriedFlux);
    end    
end

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

protPool = -ecModel.lb(strcmp(ecModel.rxns,'prot_pool_exchange'));

for i=1:numel(topEnzyme)
    [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel,topEnzyme{i});
    % See if all reactions carried flux
    [~,rIdx] = ismember(rxns,ecModel.rxns);
    carriedFlux = usageData.fluxes(rIdx) > 1e-7;
    if isscalar(find(carriedFlux))
        topUsage.protID(end+1,1)      = topEnzyme(i);
        topUsage.geneID(end+1,1)      = geneIDs(i);
        topUsage.absUsage(end+1,1)    = usageData.absUsage(topUse(i));
        topUsage.percUsage(end+1,1)   = topUsage.absUsage(end,1)/protPool*100;
        topUsage.kcat(end+1,1)        = kcat(carriedFlux);
        topUsage.source(end+1,1)      = ecModel.ec.source(idx(carriedFlux));
        topUsage.rxnID(end+1,1)       = rxns(carriedFlux);
        topUsage.rxnNames(end+1,1)    = rxnNames(carriedFlux);
        topUsage.grRules(end+1,1)     = grRules(carriedFlux);
    else
        % Add one entry for combined usage
        topUsage.protID(end+1,1)      = topEnzyme(i);
        topUsage.geneID(end+1,1)      = geneIDs(i);
        topUsage.absUsage(end+1,1)    = usageData.absUsage(topUse(i));
        topUsage.percUsage(end+1,1)   = topUsage.absUsage(end,1)/protPool*100;
        topUsage.kcat(end+1,1)        = nan;
        topUsage.source{end+1,1}      = '===';
        topUsage.rxnID{end+1,1}       = '===';
        topUsage.rxnNames{end+1,1}    = 'involved in multiple rxns, usage combined, individual rxns below';
        topUsage.grRules{end+1,1}     = '===';
        % Recalculate reaction-specific usage
        rIdx = rIdx(carriedFlux);
        enzFlux = usageData.fluxes(rIdx);
        enzMet = strcat('prot_',topEnzyme{i});
        [~, enzIdx] = ismember(enzMet,ecModel.mets);
        indUse = full(transpose(-ecModel.S(enzIdx,rIdx)).*enzFlux);
        rxnNumber = length(rIdx);
        topUsage.protID(end+1:end+rxnNumber,1)      = topEnzyme(i);
        topUsage.geneID(end+1:end+rxnNumber,1)      = geneIDs(i);
        topUsage.absUsage(end+1:end+rxnNumber,1)    = indUse;
        topUsage.percUsage(end+1:end+rxnNumber,1)   = topUsage.absUsage(end-(rxnNumber-1):end,1)/protPool*100;
        topUsage.kcat(end+1:end+rxnNumber,1)        = kcat(carriedFlux);
        topUsage.source(end+1:end+rxnNumber,1)      = ecModel.ec.source(idx(carriedFlux));
        topUsage.rxnID(end+1:end+rxnNumber,1)       = rxns(carriedFlux);
        topUsage.rxnNames(end+1:end+rxnNumber,1)    = rxnNames(carriedFlux);
        topUsage.grRules(end+1:end+rxnNumber,1)     = grRules(carriedFlux);
    end
end
usageReport.topAbsUsage     = struct2table(topUsage);
usageReport.totalUsageFlux  = protPool;
end
