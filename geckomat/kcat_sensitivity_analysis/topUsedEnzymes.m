%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function T = topUsedEnzymes(fluxes,model,conditions,name)
%
% Function that gets an ecModel_batch and a matrix of pFBA solution vectors 
% and calculate the top ten enzyme usages for every solution (conditions)
% in a mass-wise way. The results are written in the topUsedEnzymes.txt
% file and stored in the container folder.
%
% Ivan Domenzain    Last edited. 2018-09-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function topUsedEnzymes(fluxes,model,conditions,name)
    Indexes    = [];
    varNames   = {};
    outputFile = table;
    % Find the enzyme usage reactions
    usages     = find(~cellfun(@isempty,strfind(model.rxnNames,'prot_')));
    usages     = usages(1:end-1);
    usages     = fluxes(usages,:);

    for i=1:length(usages(1,:))
        %Units conversion [mmol/gDwh] -> [g prot/gDwh]
        enzUsage{i}              = usages(:,i).*model.MWs;
        %Get the mass fraction of the proteome for every protein in the mod
        enzUsage{i}              = enzUsage{i}/sum(enzUsage{i}); 
        %Sort and keep the top ten used enzymes
        [enzUsage{i},Indexes{i}] = sort( enzUsage{i},'descend');
         enzNames{i} = model.enzymes(Indexes{i});
         enzUsage{i} = (enzUsage{i}(1:10));
         enzNames{i} = enzNames{i}(1:10);
         for j=1:10
             outputFile(j,(2*i)-1) = enzNames{i}(j); 
             outputFile(j,(2*i))   = {enzUsage{i}(j)};
             varNames              = [varNames conditions(i)];
         end
    end
    %Write the top-ten used enzyme names and their percentage usages for
    %every condition on the output file
    outputFile = truncateValues(outputFile,1);
    writetable(outputFile,['../../models/' name '/data/' name '_topUsedEnzymes.txt'])
end


