%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function T = topUsedEnzymes(fluxes,model,conditions,name,writeFile,n)
%
% Function that gets an ecModel_batch and a matrix of pFBA solution vectors
% and calculate the top ten enzyme usages for every solution (conditions)
% in a mass-wise way. The results are written in the topUsedEnzymes.txt
% file and stored in the container folder.
%
%   fluxes      cell array with pFBA solution vectors
%   model       ecGEM structure used to generate pFBA solution vectors
%   conditions  cell array of strings, with names identifying each pFBA
%               solution vector
%   name        string with model abbreviation (opt, default 'ecModel')
%   writeFile   logical, whether a file should be written at the location
%               gecko/models/'name' (opt, default true)
%   n           number of top used enzymes to be listed. (opt, default 10)
%
% Ivan Domenzain    Last edited 2018-09-25
% Eduard Kerkhoven  Last edited 2018-12-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = topUsedEnzymes(fluxes,model,conditions,name,writeFile,n)

if nargin < 6
    n = 10;
end
if nargin < 5
    writeFile = true;
end
% Find the enzyme usage reactions
usages = find(~cellfun(@isempty,strfind(model.rxnNames,'prot_')));
usages = usages(1:end-1);
usages = fluxes(usages,:);
[~, nConds] = size(usages);
outputFile = cell(n,3*nConds);

for i=1:length(usages(1,:))
    %Units conversion [mmol/gDwh] -> [g prot/gDwh]
    enzUsage = usages(:,i).*model.MWs;
    %Get the mass fraction of the proteome for every protein in the mod
    enzUsage = enzUsage/sum(enzUsage);
    %Sort and keep the top ten used enzymes
    [enzUsage,Indexes] = sort(enzUsage,'descend');
    enzNames = model.enzymes(Indexes);
    enzUsage = enzUsage(1:n);
    enzNames = enzNames(1:n);
    outputFile(:,(3*i)-2)   = enzNames;
    outputFile(:,(3*i)-1)   = num2cell(enzUsage);
    outputFile(:,(3*i))     = conditions(i);
end
%Write the top-ten used enzyme names and their percentage usages for
%every condition on the output file
i=1:length(usages(1,:));
i=(3*i)-1;
outputFile = truncateValues(outputFile,i);
if writeFile
    writetable(cell2table(outputFile),['../../models/' name '/' name '_topUsedEnzymes.txt'])
end
end