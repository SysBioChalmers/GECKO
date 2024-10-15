function model = getKcatAcrossIsozymes(model)
% getKcatAcrossIsozymes
%   For reactions without kcat value (0 in model.ec.kcat), isozymes are
%   found (being based on the same reaction in the conventional GEM), that
%   do have a kcat value assigned. The mean kcat value of these isozymes
%   is then used to fill in model.ec.kcat.
%
% Input:
%   model       an ecModel in full GECKO 3 format (with ecModel.ec structure),
%               not GECKO light
%
% Output:
%   model       an ecModel with kcat values assigned to isozymes in model.ec.kcat
%
% Usage: model = getKcatAcrossIsozymes(model);

if model.ec.geckoLight
    error('Provided model is a GECKO light version, this function is not relevant for such models')
end
if all(model.ec.kcat==0)
    printOrange('WARNING: No kcat values are provided in model.ec.kcat, model remains unchanged.\n');
    return
end

noKcats     = model.ec.kcat==0;
rxnIDs      = regexprep(model.ec.rxns,'_EXP_\d+','');
noKcatID    = rxnIDs(noKcats);
yesKcatID   = rxnIDs(~noKcats);
yesKcatVal  = model.ec.kcat(~noKcats);

noKcatVal   = cellfun(@(x) strcmp(x, yesKcatID), noKcatID, 'UniformOutput', false);
noKcatVal   = cell2mat(cellfun(@(x) mean(yesKcatVal(x)), noKcatVal, 'UniformOutput', false));

newKcat     = find(~isnan(noKcatVal));
newKcatIdx  = find(noKcats);
newKcatIdx  = newKcatIdx(newKcat);
newKcat     = noKcatVal(newKcat);

model.ec.kcat(newKcatIdx) = newKcat;
model.ec.source(newKcatIdx) = {'isozymes'};
end

