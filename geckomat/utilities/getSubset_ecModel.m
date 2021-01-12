function small_ecModel = getSubset_ecModel(smallGEM,big_ecModel)
%getSubset_ecModel
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
%   small_ecModel  Enzyme-constrained version of the context-specific model
% 
% usage: small_ecModel = getSubset_ecModel(smallGEM,big_ecModel)
% 

%If smallGEM is in COBRA format convert to RAVEN-type
if isfield(smallGEM,'rules')
    smallGEM = ravenCobraWrapper(smallGEM);
end

%Check if original big_ecModel contains context-dependent protein
%constraints
measProts = ismember(big_ecModel.rxns, strcat('prot_', big_ecModel.enzymes ,'_exchange'));
if any(measProts)
    warning('The provided general ecModel contains reactions constrained by context-specific proteomics data (prot_XXXX_exchange), check if these constraints are compatible with those applied to smallGEM')
end

%Open all exchanges, this prevents removal of well-conected reactions that
%are blocked due to the imposed constraints on big_ecModel
[~,idxs]     = getExchangeRxns(big_ecModel);
if isfield(big_ecModel, 'annotation') && isfield(big_ecModel.annotation, 'defaultUB')
    big_ecModel.ub(idxs) = big_ecModel.annotation.defaultUB;
else
    big_ecModel.ub(idxs) = 1000;
end
big_ecModel.lb(idxs) = 0;

%Identify genes that are not present in smallGEM and remove all other model
%components associated to them
[~,toRemove]  = setdiff(big_ecModel.genes,smallGEM.genes);
small_ecModel = removeGenes(big_ecModel,toRemove,true,true,true);

%Remove reactions in excess (all remaining reactions in small_ecModel
%that are not represented in the context-specific GEM, mostly non 
%gene-associated reactions).
%The ecGEM reaction IDs are trimmed to remove prefix "arm_" and suffixes
%"_REV" or "No#" (where # is an integer(s)). This yields the reaction IDs
%of the original model to which they are associated.
originalRxns = smallGEM.rxns;
ecRxns_trim  = regexprep(small_ecModel.rxns, '^arm_|No\d+$|_REV($|No\d+$)', '');
toKeep       = find(ismember(ecRxns_trim, originalRxns));

%Keep enzyme-related reactions
idxs      = find(startsWith(small_ecModel.rxns,'draw_prot_') | ...
                 ismember(small_ecModel.rxns, strcat('prot_', small_ecModel.enzymes ,'_exchange')) | ...
                 strcmpi(small_ecModel.rxns,'prot_pool_exchange'));
toKeep        = [toKeep;idxs];
toRemove      = setdiff((1:numel(small_ecModel.rxns))',toKeep);
small_ecModel = removeReactions(small_ecModel,toRemove,true,true);

%obtain indexes of the enzymes that remain as pseudometabolites  in the
%reduced network
enzymes = big_ecModel.enzymes;
idxs    = find(ismember(strcat('prot_', enzymes), small_ecModel.mets));

%Correct enzyme related fields in order to remove enzymes that were removed
%from the stoichiometric matrix in the removeGenes step
f = {'enzymes', 'enzNames', 'enzGenes', 'MWs', 'sequences', 'pathways', 'concs'};
f = intersect(fieldnames(small_ecModel), f);
for i = 1:numel(f)
    small_ecModel.(f{i}) = small_ecModel.(f{i})(idxs);
end
end
 