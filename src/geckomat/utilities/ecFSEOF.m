function fseof = ecFSEOF(model,targetRxn,csRxn,nSteps,outputFile,filePath,filterG,modelAdapter)
% ecFSEOF
%   Function that runs Flux-Scanning with Enforced Objective Function (FSEOF)
%   for a specified production target.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure).
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   csRxn           rxn ID for the main carbon source uptake reaction.
%   nSteps          number of steps for suboptimal objective in FSEOF.
%                   (Optional, default 16)
%   outputFile      bolean option to save results in a file. (Optional,
%                   default false)
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
%   fseof        structure with all results. Contains the following fields:
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
%   fseof = ecFSEOF(model,targetRxn,csRxn,nSteps,outputFile,filePath,filterG,modelAdapter)

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

if nargin < 5 || isempty(outputFile)
    outputFile = false;
end

if nargin < 4 || isempty(nSteps)
    nSteps = 16;
end

% Get relevant rxn indexes
targetRxnIdx = getIndexes(model, targetRxn,'rxns');
csRxnIdx     = getIndexes(model, csRxn,'rxns');
bioRxnIdx    = getIndexes(model, params.bioRxn,'rxns');

% Standardize grRules and rxnGeneMat in model
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules        = grRules;
model.rxnGeneMat     = rxnGeneMat;

% Check carbon source uptake rate and LB defined
model = setParam(model, 'obj', params.bioRxn, 1);
sol   = solveLP(model);
if model.lb(csRxnIdx) < sol.x(csRxnIdx)
    printOrange(['WARNING: Carbon source lower bound was set to ' num2str(model.lb(csRxnIdx)) ...
        ', but the uptake rate after model optimization is ' num2str(sol.x(csRxnIdx)) '.\n'])
end

% run FSEOF analysis

% Find out the initial target production.
iniTarget = sol.x(targetRxnIdx);

% Find out the maximum theoretical yield of target reaction.
model     = setParam(model,'obj',targetRxn,1);
sol       = solveLP(model,1);
% Set to 90%% based on https://doi.org/10.1128/AEM.00115-10
maxTarget = sol.x(targetRxnIdx) * 0.9;

% Define alpha vector for suboptimal enforced objective values between
% minimal production and 90%% of the maximum theoretical yield, initialize
% fluxes, K_scores matrices and other variables
alpha     = iniTarget:((maxTarget-iniTarget)/(nSteps-1)):maxTarget;
v_matrix  = zeros(length(model.rxns),length(alpha));
% k_matrix  = zeros(length(model.rxns),length(alpha));
tolerance = 1e-4;

% % Get WT flux distribution
% model.lb(bioRxnIdx) = sol.x(bioRxnIdx) * 0.99;
% model = setParam(model, 'obj', 'prot_pool_exchange', 1);
% solWT = solveLP(model,1);
% Enforce objective flux iteratively
progressbar('Flux Scanning with Enforced Objective Function')
for i = 1:nSteps
    % Enforce the objetive flux of product formation
    model = setParam(model,'eq',targetRxnIdx,alpha(i));

    % Restore minimum biomass lb to zero and set it as objective
    model.lb(bioRxnIdx) = 0;
    model = setParam(model,'obj',params.bioRxn,1);
    sol = solveLP(model,1);

    % Relax target lb, fix biomass and minimize protein usage
    % model.lb(targetRxnIdx) = model.lb(targetRxnIdx) * 0.99;
    % model.lb(bioRxnIdx)    = sol.x(bioRxnIdx) * 0.999;
    % model = setParam(model, 'obj', 'prot_pool_exchange', 1);
    % sol = solveLP(model,1);

    % Store flux distribution
    v_matrix(:,i) = sol.x;

    % k_matrix(:,i) = v_matrix(:,i)./v_matrix(:,1);

    progressbar(i/nSteps)
end
progressbar(1) % Make sure it closes

% Take out rxns with no grRule and standard gene
withGR   = ~cellfun(@isempty,model.grRules);
stdIdx   = contains(model.grRules,'standard');
withGR(stdIdx) = 0;
rxns     = model.rxns(withGR);
v_matrix = v_matrix(withGR,:);
% k_matrix = k_matrix(withGR,:);
rxnGeneM = model.rxnGeneMat(withGR,:);

% Filter out rxns that are always zero
zero_flux = ~all(abs(v_matrix(:,1:nSteps))<=1e-2,2);
% zero_flux = sum(~isnan(k_matrix),2) > 0;
rxns      = rxns(zero_flux,:);
v_matrix  = v_matrix(zero_flux,:);
% k_matrix  = k_matrix(zero_flux,:);
rxnGeneM  = rxnGeneM(zero_flux,:);

% Identify those rxns that always increase or decrease, and only keep them
k_rxns = zeros(size(rxns));
target_rxns = logical(size(rxns));
target_type = cell(size(rxns));
for i = 1:length(rxns)
    if issorted(abs(v_matrix(i,1:nSteps)),'strictascend')
        target_rxns(i) = true;
        k_rxns(i) = abs(v_matrix(i,nSteps-1)-v_matrix(i,1))/abs(maxTarget-maxTarget/nSteps-1);
        % k_rxns(i) = abs(iniTarget-maxTarget/nSteps)/abs(v_matrix(i,nSteps)-v_matrix(i,1));
        target_type(i) = {'OE'};
    elseif issorted(abs(v_matrix(i,1:nSteps)),'strictdescend')
        target_rxns(i) = true;
        k_rxns(i) = abs(v_matrix(i,nSteps-1)-v_matrix(i,1))/abs(maxTarget-maxTarget/nSteps-1);
        if v_matrix(i,nSteps) == 0
            target_type(i) = {'KO'};
        else
            target_type(i) = {'KD'};
        end
    end
end

rxns     = rxns(target_rxns);
v_matrix = v_matrix(target_rxns,:);
% k_matrix = k_matrix(target_rxns,:);
rxnGeneM = rxnGeneM(target_rxns,:);
k_rxns   = k_rxns(target_rxns);
target_type = target_type(target_rxns);

% Order from highest to lowest
[~,order] = sort(k_rxns,'descend');
rxns      = rxns(order);
v_matrix  = v_matrix(order,:);
rxnGeneM  = rxnGeneM(order,:);
k_rxns    = k_rxns(order);
target_type = target_type(order);

ans = 1;
% % Replace remaining NaNs with 1s:
% k_matrix(isnan(k_matrix)) = 1;
% 
% % Replace any Inf value with 1000 (maximum value is ~700):
% k_matrix(k_matrix > 1000)  = 1000;
% % k_matrix(k_matrix < -1000) = 2;

% % Filter out values that are inconsistent at different alphas:
% always_down = sum(v_matrix <= 1,2) == length(alpha);
% always_up   = sum(v_matrix >= 1,2) == length(alpha);
% 
% % Identify those reactions with mixed patterns
% incons_rxns  = always_down + always_up == 0;
% 
% % Identify genes that are linked to "inconsistent rxns"
% incons_genes = sum(rxnGeneM(incons_rxns,:),1) > 0;
% 
% % Finally, inconsistent reactions are those that are not conected
% % to "inconsistent genes" from the original "inconsistent rxns" set
% incons_rxns  = sum(rxnGeneM(:,incons_genes),2) > 0;
% 
% % Keep results for the consistent rxns exclusively
% rxns     = rxns(~incons_rxns,:);
% v_matrix = v_matrix(~incons_rxns,:);
% k_matrix = k_matrix(~incons_rxns,:);
% rxnGeneM = rxnGeneM(~incons_rxns,:);

% % Get median k-score across steps and order from highest to lowest
% k_rxns    = mean(k_matrix,2);
% [~,order] = sort(k_rxns,'descend');
% rxns      = rxns(order,:);
% v_matrix  = v_matrix(order,:);
% k_matrix  = k_matrix(order,:);
% rxnGeneM  = rxnGeneM(order,:);
% k_rxns    = k_rxns(order,:);

% Create list of remaining genes and filter out any inconsistent score:
% Just those genes that are connected to the remaining rxns are
genes      = model.genes(sum(rxnGeneM,1) > 0);
k_genes    = zeros(size(genes));
cons_genes = false(size(genes));
rxnGeneM   = rxnGeneM(:,sum(rxnGeneM,1) > 0);
target_type_genes = cell(size(genes));

for i = 1:length(genes)
    % Extract all the K_scores (from rxns across alphas) conected to
    % each remaining gene
    k_set         = k_rxns(rxnGeneM(:,i) > 0); disp(i);
    % % % Check the kind of control that gene i-th exerts over its reactions
    % always_down   = sum(k_set <= (1-tolerance)) == length(k_set);
    % always_up     = sum(k_set >= (1+tolerance)) == length(k_set);
    % % Evaluate if gene is always exerting either a positive or negative
    % % control
    % cons_genes(i) = always_down + always_up == 1;
    k_genes(i)    = mean(k_set);
    target_type_genes(i) = join(unique(target_type(rxnGeneM(:,i) > 0)),', ');
end

% Keep "consistent genes"
% genes    = genes(cons_genes);
% k_genes  = k_genes(cons_genes);
% rxnGeneM = rxnGeneM(:,cons_genes);

% if filterG
%     % Filter any value between mean(alpha) and 1:
%     unchanged = (k_genes >= mean(alpha) - tolerance) + (k_genes <= 1 + tolerance) == 2;
%     genes     = genes(~unchanged);
%     k_genes   = k_genes(~unchanged);
%     rxnGeneM  = rxnGeneM(:,~unchanged);
%     % Update results for rxns-related fields (remove remaining reactions
%     % without any associated gene in rxnGeneM)
%     toKeep   = (sum(rxnGeneM,2) > 0);
%     rxns     = rxns(toKeep,:);
%     v_matrix = v_matrix(toKeep,:);
%     k_matrix = k_matrix(toKeep,:);
%     k_rxns   = k_rxns(toKeep,:);
% end

% Order genes from highest to lowest k:
[~,order] = sort(k_genes,'descend');
genes     = genes(order,:);
k_genes   = k_genes(order,:);
target_type_genes = target_type_genes(order,:);

% Create output (exclude enzyme usage reactions):
toKeep                 = ~startsWith(rxns,'usage_prot_');
rxnIdx                 = getIndexes(model,rxns(toKeep),'rxns');
geneIdx                = cellfun(@(x) find(strcmpi(model.genes,x)),genes);
fseof.flux_WT          = solWT.x;
fseof.alpha            = alpha;
fseof.v_matrix         = v_matrix(toKeep,:);
fseof.k_matrix         = k_matrix(toKeep,:);
fseof.rxnsTable(:,1)   = model.rxns(rxnIdx);
fseof.rxnsTable(:,2)   = model.rxnNames(rxnIdx);
fseof.rxnsTable(:,3)   = num2cell(k_rxns(toKeep));
fseof.rxnsTable(:,4)   = model.grRules(rxnIdx);
fseof.rxnsTable(:,5)   = constructEquations(model,rxnIdx);
fseof.geneTable(:,1)   = genes;
fseof.geneTable(:,2)   = model.geneShortNames(geneIdx);
fseof.geneTable(:,3)   = num2cell(k_genes);
fseof.geneTable(:,4)   = target_type_genes;

if outputFile
    % Write file with gene targets
    writetable(cell2table(fseof.geneTable, ...
        'VariableNames', {'gene_IDs' 'gene_names' 'K_score'}), ...
        fullfile(filePath, 'ecFSEOF_genes.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    % Write file with rxn targets
    writetable(cell2table(fseof.rxnsTable, ...
        'VariableNames', {'rxn_IDs' 'rxnNames' 'K_score' 'grRules' 'rxn_formula'}), ...
        fullfile(filePath, 'ecFSEOF_rxns.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    disp(['ecFSEOF results stored at: ' newline filePath]);
end

end
