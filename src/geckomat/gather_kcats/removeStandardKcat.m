function model = removeStandardKcat(model)
% removeStandardKcat
%   Remove the "standard" pseudoenzyme and standard kcat and MW values, as
%   they were introduced by getStandardKcat. Also standard kcat values that
%   were assigned by getStandardKcat if fillZeroKcat was set to true are
%   removed from model.ec.kcat. Both the model.ec and model.S structures
%   are modified.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%                   that has standard kcat values implemented by
%                   getStandardKcat.
%
% Output:
%   model           ecModel without standard pseudoprotein and standard
%                   kcat values
%
% Usage:
%    model = removeStandardKcat(model);

% Remove standard enzyme from ec structure
stdEnzIdx = find(strcmpi(model.ec.enzymes, 'standard'));
if ~isempty(stdEnzIdx)
    model.ec.genes(stdEnzIdx)       = [];
    model.ec.enzymes(stdEnzIdx)     = [];
    model.ec.mw(stdEnzIdx)          = [];
    model.ec.sequence(stdEnzIdx)    = [];
    if isfield(model.ec,'concs')
        model.ec.concs(stdEnzIdx)   = [];
    end
    rxnEnzIdx = find(model.ec.rxnEnzMat(:,stdEnzIdx));
    model.ec.rxns(rxnEnzIdx)        = [];
    model.ec.kcat(rxnEnzIdx)        = [];
    model.ec.source(rxnEnzIdx)      = [];
    model.ec.notes(rxnEnzIdx)       = [];
    model.ec.eccodes(rxnEnzIdx)     = [];
    model.ec.rxnEnzMat(:,stdEnzIdx) = [];
    model.ec.rxnEnzMat(rxnEnzIdx,:) = [];
end

% Remove standard kcat values and reapply kcat constraints for those
% specific reactions
stdKcatIdx = find(strcmpi(model.ec.source, 'standard'));
if ~isempty(stdKcatIdx)
    model.ec.source(stdKcatIdx)     = {''};
    model.ec.kcat(stdKcatIdx)       = 0;
    model = applyKcatConstraints(model,stdKcatIdx);
end

% Remove standard protein, usage and gene from model structure
stdMetIdx = find(strcmpi(model.mets, 'prot_standard'));
if ~isempty(stdMetIdx)
    model       = removeMets(model,stdMetIdx,false,false,false,false);
end
stdProtEx = find(strcmpi(model.rxns, 'usage_prot_standard'));
if ~isempty(stdProtEx)
    model       = removeReactions(model,stdProtEx,false,false,false);
end
stdProtGene = find(strcmpi(model.genes, 'standard'));
if ~isempty(stdProtGene)
    model       = removeGenes(model,stdProtGene,false,false,false);
end
end
