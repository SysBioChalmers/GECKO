function results = ecFSEOF(ecModel,targetRxn,csRxn,alphaLims,nSteps,file_genes,file_rxns,modelAdapter)
% ecFSEOF
%
% Function that runs Flux-scanning with Enforced Objective Function
% for a specified production target.
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure).
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   csRxn           rxn ID for the main carbon source uptake reaction.
%   alphaLims       vector of Minimum and maximum biomass scalling factors for
%                   enforced objective limits (e.g. [0.5 1]). Max value: 1.
%   nSteps          number of steps for suboptimal objective in FSEOF.
%                   (Optional, default 16)
%   file_genes      file name for results output at the genes level.
%                   (Optional, default in the 'output' sub-folder taken from
%                   modelAdapter, e.g. output/ecFSEOF_genes.tsv
%   file_rxns       file name for results output at the reactions level.
%                   (Optional, default in the 'output' sub-folder taken from
%                   modelAdapter, e.g. output/ecFSEOF_rxns.tsv)
%   modelAdapter    a loaded model adapter. (Optional, will otherwise use
%                   the default model adapter)
%
% Usage:
%   results = ecFSEOF(ecModel,targetRxn,csRxn,alphaLims)
%

if nargin < 8 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 7 || isempty(file_rxns)
    file_rxns = fullfile(params.path,'output','ecFSEOF_rxns.tsv');
end

if nargin < 6 || isempty(file_genes)
    file_genes = fullfile(params.path,'output','ecFSEOF_genes.tsv');
end

if nargin < 5 || isempty(nSteps)
    nSteps = 16;
end

if numel(alphaLims) ~= 2
    error('alphaLims parameter should be a vector of two values')
end

% Define alpha vector for suboptimal enforced objective values
alphaV  = alphaLims(1):((alphaLims(2)-alphaLims(1))/(nSteps-1)):alphaLims(2);

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
results = ecFluxScanning(ecModel,targetRxn,csRxn,alphaV,1e-4,true);

% Create gene table:
results.geneTable      = cell(length(results.genes),3);
results.geneTable(:,1) = results.genes;
results.geneTable(:,2) = results.geneNames;
results.geneTable(:,3) = num2cell(results.k_genes);

% Create rxns table (exclude enzyme usage reactions):
toKeep                 = find(~startsWith(results.rxns(:,1),'usage_prot_'));
results.k_rxns         = results.k_rxns(toKeep);
results.k_matrix       = results.k_matrix(toKeep,:);
results.v_matrix       = results.v_matrix(toKeep,:);
results.rxnsTable      = cell(length(results.k_rxns),5);
results.rxnsTable(:,1) = results.rxns(toKeep,1);
results.rxnsTable(:,2) = results.rxns(toKeep,2);
results.rxnsTable(:,3) = num2cell(results.k_rxns);
results.rxnsTable(:,4) = results.rxns(toKeep,3);
results.rxnsTable(:,5) = results.rxns(toKeep,4);

writetable(cell2table(results.geneTable, ...
    'VariableNames', {'gene_IDs' 'gene_names' 'K_score'}), ...
    file_genes, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'QuoteStrings', false);

writetable(cell2table(results.rxnsTable, ...
    'VariableNames', {'rxn_IDs' 'rxnNames' 'K_score' 'grRules' 'rxn_formula'}), ...
    file_rxns, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'QuoteStrings', false);

disp(['ecFSEOF results stored at: ' newline fileparts(file_genes)]);

% Remove redundant output fields
results = rmfield(results,'k_rxns');
results = rmfield(results,'rxns');
results = rmfield(results,'genes');
results = rmfield(results,'geneNames');
results = rmfield(results,'k_genes');
end
