function model = getKcatAcrossIsoenzymes(model)
% getKcatAcrossIsoenzymes
%   For reactions without kcat value (0 in model.ec.kcat), isoenzymes are
%   found (being based on the same reaction in the conventional GEM), that
%   do have a kcat value assigned. The mean kcat value of these isoenzymes
%   is then used to fill in model.ec.kcat.
%
% Input:
%   model       an ecModel in GECKO 3 version, not geckoLight
%
% Output:
%   model       an ecModel in GECKO 3 version with kcat values assigned to
%               isoenzymes in model.ec.kcat
%
% Usage: model = getKcatAcrossIsoenzymes(model);

if model.ec.geckoLight
    error('Provided model is a GECKO light version, this function is not relevant for such models')
end
if all(model.ec.kcat==0)
    warning('No kcat values are provided in model.ec.kcat, model remains unchanged.')
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

