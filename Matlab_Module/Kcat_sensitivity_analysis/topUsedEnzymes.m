function T = topUsedEnzymes(fluxes,model2,conditions)
MWs        = model2.MWs;
enzymes    = model2.enzymes;
enz_usages = [];
Indexes    = [];
topOne     = [];  
usages     = find(~cellfun(@isempty,strfind(model2.rxnNames,'prot_')));
usages     = usages(1:end-1);
usages     = fluxes(usages,:);
varNames   = [];
T = table;
for i=1:length(usages(1,:))
    enzUsage{i}              = usages(:,i).*model2.MWs;
    enzUsage{i}              = enzUsage{i}/sum(enzUsage{i}); 
    [enzUsage{i},Indexes{i}] = sort( enzUsage{i},'descend');
     enzNames{i} = model2.enzymes(Indexes{i});
     enzUsage{i} = (enzUsage{i}(1:10));
     enzNames{i} = enzNames{i}(1:10);
     for j=1:10
         T(j,(2*i)-1) = enzNames{i}(j); 
         T(j,(2*i))   = {enzUsage{i}(j)};
         varNames     = [varNames conditions(i)];
     end
end
writetable(T,'topUsedEnzymes.txt')
end


