function small_ecModel = getSubset_ecModel(smallGEM,big_ecModel)
%getsmall_ecModel
%  
% Generate a context-specific ecModel (strain/cell-line/tissue) by mapping 
% all components in an original context-specific GEM to a general ecModel
% (organism or species level).
% 
%   big_ecModel    Enzyme-constrained version of  general model of
%                  metabolism for a given species
%   smallGEM       Reduced model (subset of the general model) for a
%                  specific strain (microbes) or cell-line/tissues (mammals)
% 
%   ecModel        Enzyme-constrained version of the context-specific model
% 
% usage: small_ecModel = getsmall_ecModel(smallGEM,big_ecModel)
% 

%If smallGEM is in COBRA format convert to RAVEN-type
if isfield(smallGEM,'rules')
    smallGEM = ravenCobraWrapper(smallGEM);
end

%Open all exchanges, this prevents removal of well-conected reactions that
%are blocked due to the imposed constraints on big_ecModel
[~,idxs]     = getExchangeRxns(big_ecModel);
big_ecModel.ub(idxs) = 1000;
big_ecModel.lb(idxs) = 0;

%Identify genes that are not present in smallGEM and remove all other model
%components associated to them
[~,toRemove]  = setdiff(big_ecModel.genes,smallGEM.genes);
small_ecModel = removeGenes(big_ecModel,toRemove,true,true,true);

%Remove reactions in excess (all remaining reactions in small_ecModel
%that are not represented in the context-specific GEM, mostly non gene-associated reactions)
originalRxns = smallGEM.rxns;
toKeep       = [];
for rxn = originalRxns
    idxs   = find(contains(small_ecModel.rxns,rxn));
    toKeep = [toKeep;idxs];
end
%Keep enzyme-related reactions
toKeep    = unique(toKeep);
idxs      = find(contains(small_ecModel.rxns,'draw_prot_') | strcmpi(small_ecModel.rxns,'prot_pool_exchange'));
toKeep    = [toKeep;idxs];
toRemove  = setdiff((1:numel(small_ecModel.rxns))',toKeep);
small_ecModel = removeReactions(small_ecModel,toRemove,true,true);

%obtain indexes of the enzymes that remain as pseudometabolites  in the
%reduced network
enzymes = big_ecModel.enzymes;
idxs    = cellfun(@(x) find(contains(small_ecModel.mets, x)), enzymes,'UniformOutput',false);
idxs    = find(~cellfun(@isempty,idxs));

%Correct enzyme related fields in order to remove enzymes that were removed
%from the stoichiometric matrix in the removeGenes step
small_ecModel.enzymes   = small_ecModel.enzymes(idxs);
small_ecModel.enzNames  = small_ecModel.enzNames(idxs);
small_ecModel.enzGenes  = small_ecModel.enzGenes(idxs);
small_ecModel.MWs       = small_ecModel.MWs(idxs);
small_ecModel.sequences = small_ecModel.sequences(idxs);
small_ecModel.pathways  = small_ecModel.pathways(idxs);
small_ecModel.concs     = small_ecModel.concs(idxs);
end
 