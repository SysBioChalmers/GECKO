function FC = ecFSEOF(ecModel,targetRxn,csRxn,alphaLims,nSteps,filePath,filterG,modelAdapter)
% ecFSEOF
%   Function that runs Flux-Scanning with Enforced Objective Function (FSEOF)
%   for a specified production target.
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure).
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   csRxn           rxn ID for the main carbon source uptake reaction.
%   alphaLims       vector of Minimum and maximum biomass scalling factors for
%                   enforced objective limits (e.g. [0.5 0.75]). Max value: 1.
%                   Max value recomended 0.9 (Optional, default [0.5 0.75])
%   nSteps          number of steps for suboptimal objective in FSEOF.
%                   (Optional, default 16)
%   filePath        file path for results output. It will store two files:
%                   - at the genes level, ecFSEOF_genes.tsv
%                   - at the reactions level, ecFSEOF_rxns.tsv
%                   (Optional, default in the 'output' sub-folder taken from
%                   modelAdapter, e.g. output/ecFSEOF_rxns.tsv)
%   filterG         logical value. TRUE if genes K_scores results should be
%                   filtered according to the alpha vector distribution
%                   (Optional, default true)
%   modelAdapter    a loaded model adapter. (Optional, will otherwise use
%                   the default model adapter)
%
% Output:
%   FC           structure with all results. Contains the following fields:
%                - flux_WT:   flux distribution for the WT strain (100% of
%                             carbon towards growth).
%                - alpha:     scalling factors used for enforced objetive
%                             limits
%                - v_matrix:  fluxes for each reaction in rxnsTable.rxns
%                             and each alpha.
%                - k_matrix:  fold changes for each reaction in
%                             rxnsTable.rxns and each alpha.
%                - rxnsTable: a list with all reactions with fluxes that
%                             change consistently as target production 
%                             increases.
%                             Contains: ID, name, k_score, gene rule, equation
%                - geneTable: a list with all selected targets that
%                             increase production.
%                             Contains: gene, shortName, k_score
%
% Usage:
%   FC = ecFSEOF(ecModel,targetRxn,csRxn,alphaLims,nSteps,filePath,filterG,modelAdapter)

if nargin < 8 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 7 || isempty(filterG)
    filterG = true;
end

if nargin < 6 || isempty(filePath)
    filePath = fullfile(params.path,'output');
end

if nargin < 5 || isempty(nSteps)
    nSteps = 16;
end

if nargin < 4 || isempty(alphaLims)
    alphaLims = [0.5 0.75];
end

if numel(alphaLims) ~= 2
    error('alphaLims parameter should be a vector of two values')
end

% Standardize grRules and rxnGeneMat in model
[grRules,rxnGeneMat] = standardizeGrRules(ecModel,true);
ecModel.grRules      = grRules;
ecModel.rxnGeneMat   = rxnGeneMat;

% Check carbon source uptake rate and LB
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel, 1);
csRxnIdx = strcmpi(ecModel.rxns,csRxn);

if sol.x(csRxnIdx) < ecModel.lb(csRxnIdx)
    printOrange('WARNING: Carbon source LB and uptake rate are not equal.')
end

% run FSEOF analysis

% Define alpha vector for suboptimal enforced objective values and
% initialize fluxes and K_scores matrices
alpha = alphaLims(1):((alphaLims(2)-alphaLims(1))/(nSteps-1)):alphaLims(2);
v_matrix = zeros(length(ecModel.rxns),length(alpha));
k_matrix = zeros(length(ecModel.rxns),length(alpha));
enzRxns  = contains(ecModel.rxns,'usage_prot_');
tolerance = 1e-4;

% Simulate WT (100% growth) and forced (X% growth and the rest towards product):
flux_WT = getFluxTarget(ecModel,params.bioRxn,targetRxn,1);

% Set values which are under solver tolerance
flux_WT(flux_WT < 0 & ~enzRxns) = 0;
flux_WT(flux_WT > 0 & enzRxns) = 0;

progressbar('Flux Scanning with Enforced Objective Function')
for i = 1:length(alpha)
    flux_MAX = getFluxTarget(ecModel,params.bioRxn,targetRxn,alpha(i));
    % Set values which are under solver tolerance
    flux_MAX(flux_MAX < 0 & ~enzRxns) = 0;
    flux_MAX(flux_MAX > 0 & enzRxns) = 0;
    v_matrix(:,i) = flux_MAX;
    k_matrix(:,i) = flux_MAX./flux_WT;
    progressbar(i/length(alpha))
end
progressbar(1) % Make sure it closes 

% take out rxns with no grRule and standard gene 
withGR   = ~cellfun(@isempty,ecModel.grRules);
stdPos   = contains(ecModel.grRules,'standard');
withGR(stdPos) = 0;
rxns     = ecModel.rxns(withGR);
v_matrix = v_matrix(withGR,:);
k_matrix = k_matrix(withGR,:);
rxnGeneM = ecModel.rxnGeneMat(withGR,:);

% Filter out rxns that are always zero -> k=0/0=NaN:
non_nan  = sum(~isnan(k_matrix),2) > 0;
rxns     = rxns(non_nan,:);
v_matrix = v_matrix(non_nan,:);
k_matrix = k_matrix(non_nan,:);
rxnGeneM = rxnGeneM(non_nan,:);

% Replace remaining NaNs with 1s:
k_matrix(isnan(k_matrix)) = 1;

% Replace any Inf value with 1000 (maximum value is ~700):
k_matrix(k_matrix > 1000)  = 1000;
% k_matrix(k_matrix < -1000) = 2;

% Filter out values that are inconsistent at different alphas:
always_down = sum(k_matrix <= 1,2) == length(alpha);
always_up   = sum(k_matrix >= 1,2) == length(alpha);

% Identify those reactions with mixed patterns
incons_rxns  = always_down + always_up == 0;

% Identify genes that are linked to "inconsistent rxns"
incons_genes = sum(rxnGeneM(incons_rxns,:),1) > 0;

% Finally, inconsistent reactions are those that are not conected
% to "inconsistent genes" from the original "inconsistent rxns" set
incons_rxns  = sum(rxnGeneM(:,incons_genes),2) > 0;

% Keep results for the consistent rxns exclusively
rxns     = rxns(~incons_rxns,:);
v_matrix = v_matrix(~incons_rxns,:);
k_matrix = k_matrix(~incons_rxns,:);
rxnGeneM = rxnGeneM(~incons_rxns,:);

% Get median k-score across steps and order from highest to lowest
k_rxns    = mean(k_matrix,2);
[~,order] = sort(k_rxns,'descend');
rxns      = rxns(order,:);
v_matrix  = v_matrix(order,:);
k_matrix  = k_matrix(order,:);
rxnGeneM  = rxnGeneM(order,:);
k_rxns    = k_rxns(order,:);

% Create list of remaining genes and filter out any inconsistent score:
% Just those genes that are connected to the remaining rxns are
genes      = ecModel.genes(sum(rxnGeneM,1) > 0);
k_genes    = zeros(size(genes));
cons_genes = false(size(genes));
rxnGeneM   = rxnGeneM(:,sum(rxnGeneM,1) > 0);

for i = 1:length(genes)
    % Extract all the K_scores (from rxns across alphas) conected to
    % each remaining gene
    k_set         = k_rxns(rxnGeneM(:,i) > 0);
    % Check the kind of control that gene i-th exerts over its reactions
    always_down   = sum(k_set <= (1-tolerance)) == length(k_set);
    always_up     = sum(k_set >= (1+tolerance)) == length(k_set);
    % Evaluate if gene is always exerting either a positive or negative
    % control
    cons_genes(i) = always_down + always_up == 1;
    k_genes(i)    = mean(k_set);
end

% Keep "consistent genes"
genes    = genes(cons_genes);
k_genes  = k_genes(cons_genes);
rxnGeneM = rxnGeneM(:,cons_genes);

if filterG
    % Filter any value between mean(alpha) and 1:
    unchanged = (k_genes >= mean(alpha) - tolerance) + (k_genes <= 1 + tolerance) == 2;
    genes     = genes(~unchanged);
    k_genes   = k_genes(~unchanged);
    rxnGeneM  = rxnGeneM(:,~unchanged);
    % Update results for rxns-related fields (remove remaining reactions
    % without any associated gene in rxnGeneM)
    toKeep   = (sum(rxnGeneM,2) > 0);
    rxns     = rxns(toKeep,:);
    v_matrix = v_matrix(toKeep,:);
    k_matrix = k_matrix(toKeep,:);
    k_rxns   = k_rxns(toKeep,:);
end

% Order genes from highest to lowest k:
[~,order] = sort(k_genes,'descend');
genes     = genes(order,:);
k_genes   = k_genes(order,:);

% Create output (exclude enzyme usage reactions):
toKeep                 = ~startsWith(rxns,'usage_prot_');
rxnIdx                 = getIndexes(ecModel,rxns(toKeep),'rxns');
geneIdx                = cellfun(@(x) find(strcmpi(ecModel.genes,x)),genes);
FC.flux_WT             = flux_WT;
FC.alpha               = alpha;
FC.v_matrix            = v_matrix(toKeep,:);
FC.k_matrix            = k_matrix(toKeep,:);
FC.rxnsTable(:,1)      = ecModel.rxns(rxnIdx);
FC.rxnsTable(:,2)      = ecModel.rxnNames(rxnIdx);
FC.rxnsTable(:,3)      = num2cell(k_rxns(toKeep));
FC.rxnsTable(:,4)      = ecModel.grRules(rxnIdx);
FC.rxnsTable(:,5)      = constructEquations(ecModel,rxnIdx);
FC.geneTable(:,1)      = genes;
FC.geneTable(:,2)      = ecModel.geneShortNames(geneIdx);
FC.geneTable(:,3)      = num2cell(k_genes);

writetable(cell2table(FC.geneTable, ...
    'VariableNames', {'gene_IDs' 'gene_names' 'K_score'}), ...
    fullfile(filePath, 'ecFSEOF_genes.tsv'), ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'QuoteStrings', false);

writetable(cell2table(FC.rxnsTable, ...
    'VariableNames', {'rxn_IDs' 'rxnNames' 'K_score' 'grRules' 'rxn_formula'}), ...
    fullfile(filePath, 'ecFSEOF_rxns.tsv'), ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'QuoteStrings', false);

disp(['ecFSEOF results stored at: ' newline filePath]);
end
