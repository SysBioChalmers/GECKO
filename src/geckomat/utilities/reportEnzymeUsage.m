function usageReport = reportEnzymeUsage(ecModel, usageData, varargin)
% reportEnzymeUsage  Summarize the results from enzymeUsage.
%
% Ranks the enzymes from enzymeUsage into two tables: those with the
% highest capacity usage, and those with the highest absolute usage.
%
% Parameters
% ----------
% ecModel : struct
%     a GECKO3 ecModel.
% usageData : struct
%     output from enzymeUsage.
%
% Name-Value Arguments
% --------------------
% highCapUsage : double
%     minimum fraction of enzyme capacity usage to be listed as high usage
%     (default 0.9, i.e. at least 90% capacity usage).
% topAbsUsage : double
%     how many top enzymes by absolute usage to list (default 10). Use 0
%     or Inf to list all enzymes instead. Enzymes whose reactions carry no
%     flux are never listed here, as they do not contribute to the
%     solution.
%
% Returns
% -------
% usageReport : struct
%     highCapUsage : table
%         enzymes at or above the highCapUsage threshold.
%     topAbsUsage : table
%         the topAbsUsage enzymes with the highest absolute usage.
%     totalUsageFlux : double
%         total protein pool used, i.e. the upper bound of prot_pool_exchange.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     usageReport = reportEnzymeUsage(ecModel, usageData, highCapUsage, topAbsUsage);
%     usageReport = reportEnzymeUsage(ecModel, usageData, 'topAbsUsage', 20);
%
% See also
% --------
% enzymeUsage

p = parseGECKOargs(varargin, { ...
    'highCapUsage', []; ...
    'topAbsUsage',  []});
highCapUsage = p.highCapUsage;
topAbsUsage  = p.topAbsUsage;
if isempty(highCapUsage); highCapUsage = 0.9; end
if isempty(topAbsUsage);  topAbsUsage  = 10;  end
% 0 or Inf means "all enzymes", per the docstring above -- and any topAbsUsage at or
% above the enzyme count means the same thing in practice, but unlike those two special
% values it used to reach the indexing below unclamped: usageData.protID(topUse(1:N))
% with N > numel(topUse) is an out-of-bounds index, so any model with fewer than 10
% enzymes (ecTestGEM has 5) crashed on the default call.
if topAbsUsage == 0 || isinf(topAbsUsage) || topAbsUsage > numel(usageData.protID)
    topAbsUsage = numel(usageData.protID);
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

protPool = ecModel.ub(strcmp(ecModel.rxns,'prot_pool_exchange'));

for i=1:numel(topEnzyme)
    [rxns, kcat, idx, rxnNames, grRules] = getReactionsFromEnzyme(ecModel,topEnzyme{i});
    % See if all reactions carried flux
    [~,rIdx] = ismember(rxns,ecModel.rxns);
    carriedFlux = usageData.fluxes(rIdx) > 1e-7;
    if ~any(carriedFlux)
        % None of this enzyme's reactions carry flux -- skip it rather than
        % padding the table with a placeholder row for an enzyme that is not
        % actually active.
        continue
    elseif isscalar(find(carriedFlux))
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
