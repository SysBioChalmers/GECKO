function smallEcModel = getSubsetEcModel(bigEcModel,smallGEM)
% getSubset_ecModel
%   Generate a context-specific ecModel (strain/cell-line/tissue) by
%   mapping reactions and genes from a conventional context-specific GEM to
%   a general ecModel (bigEcModel), to yield a context-specific ecModel.
% 
% Input:
%   bigEcModel  enzyme-constrained version of conventional model of
%                  metabolism for a given species
%   smallGEM       Reduced model (subset of the general model) for a
%                  specific strain (microbes) or cell-line/tissues (mammals)
% 
% Output:
%   smallEcModel  Enzyme-constrained version of the context-specific model
% 
% Usage: smallEcModel = getSubsetEcModel(smallGEM,bigEcModel)

% Check if original bigEcModel contains context-dependent protein constraints
if any(bigEcModel.lb(startsWith(bigEcModel.rxns,'usage_prot_')) ~= -1000)
    disp(['The bigEcModel is constraint by protein concentrations that are ' ...
             'likely not relevant in the constructed smallEcModel.'])
end

% Remove genes (and associated reactions) that are absent in smallGEM
genesToRemove  = setdiff(bigEcModel.genes,smallGEM.genes);
smallEcModel = removeGenes(bigEcModel,genesToRemove,true,true,false);

% Remove genes from ec-structure
enzToRemove = ismember(smallEcModel.ec.genes,genesToRemove);
smallEcModel.ec.genes(enzToRemove)       = [];
smallEcModel.ec.enzymes(enzToRemove)     = [];
smallEcModel.ec.concs(enzToRemove)       = [];
smallEcModel.ec.mw(enzToRemove)          = [];
smallEcModel.ec.sequence(enzToRemove)    = [];
smallEcModel.ec.rxnEnzMat(:,enzToRemove) = [];

% Remove any reaction (except prot_ associated) that is absent from smallGEM
trimRxns = regexprep(smallEcModel.rxns,'_REV|_EXP_\d+','');
protRxns = startsWith(smallEcModel.rxns,{'usage_prot_','prot_pool_exchange'});
keepRxns = ismember(trimRxns,smallGEM.rxns) | protRxns;
smallEcModel = removeReactions(smallEcModel,~keepRxns, true, true, true);

% Remove reactions from ec-structure
ecRxnsToRemove = ~ismember(smallEcModel.ec.rxns,smallEcModel.rxns);

smallEcModel.ec.rxns(ecRxnsToRemove)        = [];
smallEcModel.ec.kcat(ecRxnsToRemove)        = [];
smallEcModel.ec.source(ecRxnsToRemove)      = [];
smallEcModel.ec.notes(ecRxnsToRemove)       = [];
smallEcModel.ec.eccodes(ecRxnsToRemove)     = [];
smallEcModel.ec.rxnEnzMat(ecRxnsToRemove,:) = [];

smallEcModel = removeReactions(smallEcModel,~keepRxns, true, true, true);
end
 