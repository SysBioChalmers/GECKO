function FC = ecFlux_scanning(ecModel,targetRxn,csRxn,alpha,tolerance,filterG)
%ecFlux_scanning
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure).
%   targetRxn       rxn ID for the production target reaction, a exchange
%                   reaction is recommended.
%   csRxn           rxn ID for the main carbon source uptake reaction.
%   alpha           scalling factor for production yield for enforced objective
%                   limits
%   tolerance       numerical tolerance for fixing bounds.
%                   (Optional, defaul 1E-4)
%   filterG         logical value. TRUE if genes K_scores results should be
%                   filtered according to the alpha vector distribution
%                   (Optional, defaul false)
%
% Usage:
%   FC = ecFlux_scanning(model,targetRxn,csRxn,alpha,tolerance,filterG)

if nargin < 6 || isempty(filterG)
    filterG = false;
end

if nargin < 5 || isempty(tolerance)
    tolerance = 1E-4;
end

% Simulate WT (100% growth):
[~, FC.flux_WT] = getFluxTarget(ecModel,targetRxn,csRxn);

% Set to zero values which are under solver tolerance
FC.flux_WT(abs(FC.flux_WT) < 1e-8) = 0; 

% Simulate forced (X% growth and the rest towards product) based on yield:
FC.alpha = alpha;

% Initialize fluxes and K_scores matrices
v_matrix = zeros(length(ecModel.rxns),length(alpha));
k_matrix = zeros(length(ecModel.rxns),length(alpha));
for i = 1:length(alpha)
    %disp(['Iteration #' num2str(i)])
    [~, FC.flux_MAX] = getFluxTarget(ecModel,targetRxn,csRxn,alpha(i));
    % Set to zero values which are under solver tolerance
    FC.flux_MAX(abs(FC.flux_MAX) < 1e-8) = 0;
    v_matrix(:,i) = FC.flux_MAX;
    k_matrix(:,i) = FC.flux_MAX./FC.flux_WT;
end

% Take out rxns with no grRule:
withGR   = ~cellfun(@isempty,ecModel.grRules);

% Generate rxn equations:
rxnEqs   = constructEquations(ecModel,withGR,true);
v_matrix = v_matrix(withGR,:);
k_matrix = k_matrix(withGR,:);
rxnGeneM = ecModel.rxnGeneMat(withGR,:);
FC.rxns  = [ecModel.rxns(withGR),ecModel.rxnNames(withGR),ecModel.grRules(withGR),rxnEqs];

% Filter out rxns that are always zero -> k=0/0=NaN:
non_nan  = sum(~isnan(k_matrix),2) > 0;
v_matrix = v_matrix(non_nan,:);
k_matrix = k_matrix(non_nan,:);
rxnGeneM = rxnGeneM(non_nan,:);
FC.rxns  = FC.rxns(non_nan,:);

% Replace remaining NaNs with 1s:
k_matrix(isnan(k_matrix)) = 1;

% Replace any Inf value with 1000 (maximum value is ~700):
k_matrix(k_matrix>1000) = 1000;
k_matrix(k_matrix<-1000) = 1000;

% Filter out values that are inconsistent at different alphas:
always_down  = sum(k_matrix <= 1,2) == length(alpha);
always_up    = sum(k_matrix >= 1,2) == length(alpha);

% Identify those reactions with mixed patterns
incons_rxns  = always_down + always_up == 0;

% Identify genes that are linked to "inconsistent rxns"
incons_genes = sum(rxnGeneM(incons_rxns,:),1) > 0;

% Finally, inconsistent reactions are those that are not conected
% to "inconsistent genes" from the original "inconsistent rxns" set
incons_rxns  = sum(rxnGeneM(:,incons_genes),2) > 0;

% Keep results for the consistent rxns exclusively
v_matrix     = v_matrix(~incons_rxns,:);
k_matrix     = k_matrix(~incons_rxns,:);
rxnGeneM     = rxnGeneM(~incons_rxns,:);
FC.rxns      = FC.rxns(~incons_rxns,:);

% Get median k score across steps
FC.k_rxns   = mean(k_matrix,2);

% Order from highest to lowest median k_score (across alphas)
[~,order]   = sort(FC.k_rxns,'descend');
FC.k_rxns   = FC.k_rxns(order,:);
FC.v_matrix = v_matrix(order,:);
FC.k_matrix = k_matrix(order,:);
rxnGeneM    = rxnGeneM(order,:);
FC.rxns     = FC.rxns(order,:);

% Create list of remaining genes and filter out any inconsistent score:
% Just those genes that are connected to the remaining rxns are
FC.genes     = ecModel.genes(sum(rxnGeneM,1) > 0);
FC.geneNames = ecModel.geneShortNames(sum(rxnGeneM,1) > 0);
FC.k_genes   = zeros(size(FC.genes));
cons_genes   = false(size(FC.genes));
rxnGeneM     = rxnGeneM(:,sum(rxnGeneM,1) > 0);

for i = 1:length(FC.genes)
    % Extract all the K_scores (from rxns across alphas) conected to
    % each remaining gene
    k_set         = FC.k_rxns(rxnGeneM(:,i) > 0);
    % Check the kind of control that gene i-th exerts over its reactions
    always_down   = sum(k_set <= (1-tolerance)) == length(k_set);
    always_up     = sum(k_set >= (1+tolerance)) == length(k_set);
    % Evaluate if gene is always exerting either a positive or negative
    % control
    cons_genes(i) = always_down + always_up == 1;
    FC.k_genes(i) = mean(k_set);
end
% Keep "consistent genes"
FC.genes     = FC.genes(cons_genes);
FC.geneNames = FC.geneNames(cons_genes);
FC.k_genes   = FC.k_genes(cons_genes);
rxnGeneM     = rxnGeneM(:,cons_genes);

if filterG
    % Filter any value between mean(alpha) and 1:
    unchanged    = (FC.k_genes >= mean(alpha) - tolerance) + (FC.k_genes <= 1 + tolerance) == 2;
    FC.genes     = FC.genes(~unchanged);
    FC.geneNames = FC.geneNames(~unchanged);
    FC.k_genes   = FC.k_genes(~unchanged);
    rxnGeneM     = rxnGeneM(:,~unchanged);
    % Update results for rxns-related fields (remove remaining reactions
    % without any associated gene in rxnGeneM)
    toKeep      = (sum(rxnGeneM,2) > 0);
    FC.k_rxns   = FC.k_rxns(toKeep,:);
    FC.v_matrix = v_matrix(toKeep,:);
    FC.k_matrix = k_matrix(toKeep,:);
    FC.rxns     = FC.rxns(toKeep,:);
end

% Order from highest to lowest k:
[~,order]    = sort(FC.k_genes,'descend');
FC.genes     = FC.genes(order,:);
FC.geneNames = FC.geneNames(order,:);
FC.k_genes   = FC.k_genes(order,:);

end
