function fseof = ecFSEOF(model,prodTargetRxn,csRxn,nSteps,outputFile,filePath,modelAdapter)
% ecFSEOF
%   Function that runs Flux-Scanning with Enforced Objective Function (FSEOF)
%   for a specified production target.
%
% Input:
%   model          an ecModel in GECKO 3 format (with ecModel.ec structure).
%   prodTargetRxn  rxn ID for the production target reaction, a exchange
%                  reaction is recommended.
%   csRxn          rxn ID for the main carbon source uptake reaction.
%   nSteps         number of steps for suboptimal objective in FSEOF.
%                  (Optional, default 16)
%   outputFile     bolean option to save results in a file. (Optional,
%                  default false)
%   filePath       file path for results output. It will store two files:
%                  - at the genes level, ecFSEOF_genes.tsv
%                  - at the reactions level, ecFSEOF_rxns.tsv
%                  (Optional, default in the 'output' sub-folder taken from
%                  modelAdapter, e.g. output/ecFSEOF_rxns.tsv)
%   modelAdapter   a loaded model adapter. (Optional, will otherwise use
%                  the default model adapter)
%
% Output:
%   fseof   an structure with all results. Contains the following fields:
%           - alpha:             target production used for enforced 
%                                objetive limits (from minimum to maximum 
%                                production)
%           - v_matrix:          fluxes for each target reaction predicted
%                                and each alpha.
%           - rxnTargets:        a list with all reactions with fluxes that
%                                change consistently as target production
%                                increases.
%                                Contains: ID, name, slope, gene rule, and
%                                equation
%           - transportTargets:  a list with all transport reactions with
%                                fluxes that change consistently as target
%                                production increases.
%                                Contains: ID, name, slope, gene rule, and
%                                equation
%           - geneTargets:       a list with all selected targets that
%                                increase production.
%                                Contains: gene, shortName, slope, action,
%                                and essentiality
%
% Usage:
%   fseof = ecFSEOF(model,prodTargetRxn,csRxn,nSteps,outputFile,filePath,filterG,modelAdapter)

if nargin < 7 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 5 || isempty(outputFile)
    outputFile = false;
end

if nargin < 6 || isempty(filePath)
    filePath = fullfile(params.path,'output');
end

if nargin < 4 || isempty(nSteps)
    nSteps = 16;
end

% Get relevant rxn indexes
prodTargetRxnIdx = getIndexes(model, prodTargetRxn,'rxns');
csRxnIdx         = getIndexes(model, csRxn,'rxns');
bioRxnIdx        = getIndexes(model, params.bioRxn,'rxns');

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
iniTarget = sol.x(prodTargetRxnIdx);

% Find out the maximum theoretical yield of target reaction.
model     = setParam(model,'obj',prodTargetRxn,1);
sol       = solveLP(model,1);
% Set to 90%% based on https://doi.org/10.1128/AEM.00115-10
maxTarget = sol.x(prodTargetRxnIdx) * 0.9;

% Define alpha vector for suboptimal enforced objective values between
% minimal production and 90%% of the maximum theoretical yield, initialize
% fluxes matriz
alpha    = iniTarget:((maxTarget-iniTarget)/(nSteps-1)):maxTarget;
v_matrix = zeros(length(model.rxns),length(alpha));

% Enforce objective flux iteratively
progressbar('Flux Scanning with Enforced Objective Function')
for i = 1:nSteps
    % Enforce the objetive flux of product formation
    model = setParam(model,'eq',prodTargetRxnIdx,alpha(i));

    % Restore minimum biomass lb to zero and set it as objective
    model.lb(bioRxnIdx) = 0;
    model = setParam(model,'obj',params.bioRxn,1);
    sol   = solveLP(model,1);

    % Store flux distribution
    v_matrix(:,i) = sol.x;

    progressbar(i/nSteps)
end
progressbar(1) % Make sure it closes

% Take out rxns with no grRule and standard gene
withGR         = ~cellfun(@isempty,model.grRules);
stdIdx         = contains(model.grRules,'standard');
withGR(stdIdx) = 0;
target_rxns    = model.rxns(withGR);
v_matrix       = v_matrix(withGR,:);
rxnGeneM       = model.rxnGeneMat(withGR,:);

% Filter out rxns that are always zero
zero_flux   = ~all(abs(v_matrix(:,1:nSteps))==0,2);
target_rxns = target_rxns(zero_flux,:);
v_matrix    = v_matrix(zero_flux,:);
rxnGeneM    = rxnGeneM(zero_flux,:);

% Identify those rxns that always increase or decrease, and calculate the
% slope as the difference in the flux when enforce objetive target 
% production is set to 90%% of the maximum teorethical yield 
% << v_matrix(i,nSteps-1) >> and the flux when the enforce objetive target
% production is set to the minimum << v_matrix(i,1) >> for the reaction i, 
% divided by maxTarget-maxTarget/nSteps-1.
slope_rxns  = zeros(size(target_rxns));
targets     = logical(size(target_rxns));
target_type = cell(size(target_rxns));
for i = 1:length(target_rxns)
    if issorted(abs(v_matrix(i,1:nSteps)),'strictascend')
        % Those reactions that always increase while enforcing target
        % production are suggested for Over Expression
        targets(i) = true;
        slope_rxns(i)  = abs(v_matrix(i,nSteps-1)-v_matrix(i,1))/abs(maxTarget-maxTarget/nSteps-1);
        target_type(i) = {'OE'};
    elseif issorted(abs(v_matrix(i,1:nSteps)),'strictdescend')
        % Those reactions that always decrease while enforcing target
        % production are suggested for KnockDown or KnockOut. KO are those
        % reactions which have zero flux when enforcing target production 
        % to 90%% of the maximum theoretical yield. 
        targets(i) = true;
        slope_rxns(i)  = abs(v_matrix(i,nSteps-1)-v_matrix(i,1))/abs(maxTarget-maxTarget/nSteps-1);
        if v_matrix(i,nSteps) == 0
            target_type(i) = {'KO'};
        else
            target_type(i) = {'KD'};
        end
    end
end

% Only keep those reaction that shows an increase or decrease pattern.
target_rxns = target_rxns(targets);
v_matrix    = v_matrix(targets,:);
rxnGeneM    = rxnGeneM(targets,:);
slope_rxns  = slope_rxns(targets);
target_type = target_type(targets);

% Order from highest to lowest slope
[~,order]   = sort(slope_rxns,'descend');
target_rxns = target_rxns(order);
v_matrix    = v_matrix(order,:);
rxnGeneM    = rxnGeneM(order,:);
slope_rxns  = slope_rxns(order);
target_type = target_type(order);

% Filter out reactions with slope < quantile
non_zero_slope  = slope_rxns > quantile(slope_rxns,0.75);
target_rxns     = target_rxns(non_zero_slope);
v_matrix        = v_matrix(non_zero_slope,:);
rxnGeneM        = rxnGeneM(non_zero_slope,:);
slope_rxns      = slope_rxns(non_zero_slope);
target_type     = target_type(non_zero_slope);

% Create gene list of those connected to the remaining rxns
genes             = model.genes(sum(rxnGeneM,1) > 0);
slope_genes       = zeros(size(genes));
rxnGeneM          = rxnGeneM(:,sum(rxnGeneM,1) > 0);
target_type_genes = cell(size(genes));
essentiality      = cell(size(genes));

% Validate for gene essentiality
progressbar('Checking for gene essentiality')
for i = 1:numel(genes)

    % Block protein usage to KO al the reactions associated to it.
    usage_rxn_idx = strcmpi(model.ec.genes, genes{i});
    usage_rxn = strcat('usage_prot_', model.ec.enzymes(usage_rxn_idx));
    tempModel = setParam(model, 'eq', usage_rxn, 0);
    solKO     = solveLP(tempModel);
    % Check if no feasible solution was found
    if solKO.stat == -1 || solKO.x(bioRxnIdx) < 1e-8
        essentiality(i) = {'essential'};
    end

    % Since a gene can be involved in multiple reactions, multiple
    % engineering manipulations can be suggested for the same gene
    % e.g. (OE and KD). So, reconcilie them, and report only one.
    
    % Get reaction index, in targets set, for each gene
    rxns_for_gene = find(rxnGeneM(:,i) > 0);
    actions = unique(target_type(rxns_for_gene));
    if numel(actions) > 1
        % If any OE is suggested give the highest priority over KD or KO
        if any(startsWith(actions, 'OE'))
            target_type_genes(i) = {'OE'};
            % Otherwise change to KO if not essential
        elseif ~any(startsWith(actions, 'OE')) && ~isequal(essentiality(i),{'essential'})
            target_type_genes(i) = {'KO'};
        else
            target_type_genes(i) = {'KD'};
        end
    else
        target_type_genes(i) = actions;
    end
    % Extract all the slope (from rxns across alphas) conected to
    % each remaining gene
    slope_set      = slope_rxns(rxns_for_gene);
    slope_genes(i) = mean(slope_set);

    progressbar(i/numel(genes))
end
progressbar(1) % Make sure it closes

% Order genes from highest to lowest slope:
[~,order]         = sort(slope_genes,'descend');
genes             = genes(order,:);
slope_genes       = slope_genes(order,:);
target_type_genes = target_type_genes(order,:);
essentiality      = essentiality(order,:);

% Create output:

% Report values used to enforce objective target production
fseof.alpha = alpha;

% Exclude enzyme usage reactions
toKeep = ~startsWith(target_rxns,'usage_prot_');
rxnIdx = getIndexes(model,target_rxns(toKeep),'rxns');

% Report the fluxes at each alpha
fseof.v_matrix = array2table(v_matrix(toKeep,:), ...
    'VariableNames', string(alpha), 'RowNames', model.rxns(rxnIdx));

% Report reaction targets
fseof.rxnTargets(:,1) = model.rxns(rxnIdx);
fseof.rxnTargets(:,2) = model.rxnNames(rxnIdx);
fseof.rxnTargets(:,3) = num2cell(slope_rxns(toKeep));
fseof.rxnTargets(:,4) = model.grRules(rxnIdx);
fseof.rxnTargets(:,5) = constructEquations(model,rxnIdx);

% Identify transport reactions; defined as reactions involving (at least) 
% one metabolite name in more than one compartment. 
transportRxns = false(numel(rxnIdx),1);
for i = 1:numel(rxnIdx)
    % Get the involved metabolites in each reaction
    mets = model.metNames(model.S(:,rxnIdx(i))~=0);
    % Remove enzyme metabolite
    mets(startsWith(mets,'prot_')) = [];
    % Validate if it is a transport reaciton
    transportRxns(i) = numel(mets) ~= numel(unique(mets));
end

% Create separate table for transport reactions
fseof.transportTargets = fseof.rxnTargets(transportRxns,:);
fseof.rxnTargets(transportRxns,:) = [];

% Report gene targets
geneIdx              = cellfun(@(x) find(strcmpi(model.genes,x)),genes);
fseof.geneTargets(:,1) = genes;
fseof.geneTargets(:,2) = model.geneShortNames(geneIdx);
fseof.geneTargets(:,3) = num2cell(slope_genes);
fseof.geneTargets(:,4) = target_type_genes;
fseof.geneTargets(:,5) = essentiality;

% Save results in a file if defined
if outputFile
    % Write file with gene targets
    writetable(cell2table(fseof.geneTargets, ...
        'VariableNames', {'gene_IDs' 'gene_names' 'slope' 'action' 'essentiality'}), ...
        fullfile(filePath, 'ecFSEOF_genes.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    % Write file with rxn targets
    writetable(cell2table(fseof.rxnTargets, ...
        'VariableNames', {'rxn_IDs' 'rxnNames' 'slope' 'grRules' 'rxn_formula'}), ...
        fullfile(filePath, 'ecFSEOF_rxns.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    % Write file with transport targets
    writetable(cell2table(fseof.transportTargets, ...
        'VariableNames', {'rxn_IDs' 'rxnNames' 'slope' 'grRules' 'rxn_formula'}), ...
        fullfile(filePath, 'ecFSEOF_transporter.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    disp(['ecFSEOF results stored at: ' newline filePath]);
end

end
