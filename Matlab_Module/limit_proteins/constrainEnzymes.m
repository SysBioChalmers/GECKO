%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainEnzymes(model,Ptotal,sigma,pIDs,data)
% 
% Benjam?n J. S?nchez. Last edited: 2018-03-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,enzUsages,modifications] = constrainEnzymes(model,Ptot,sigma,pIDs,data)

%Current values:
f       = 0.4461; %Yeast 7.6 (all enzymes) [g(Pmodel)/g(Ptot)]
Pbase   = 0.4005; %Value from biomass comp. (F?rster data @ 0.1 1/h)

%No UB will be changed if no data is available -> pool = all enzymes(FBAwMC)
if nargin == 3
    pIDs = {};
    data = [];
end
%Remove zeros or negative values
data = cleanDataset(data);
%Assign concentrations as UBs [mmol/gDW]:
model.concs = nan(size(model.enzymes));      %OBS: min value is zero!!
disp('Matching data to enzymes in model...')
for i = 1:length(model.enzymes)
    match = false;
    for j = 1:length(pIDs)
        if strcmpi(pIDs{j},model.enzymes{i}) && ~match
            model.concs(i) = data(j)*model.MWs(i); %g/gDW
            rxn_name       = ['prot_' model.enzymes{i} '_exchange'];
            pos            = strcmpi(rxn_name,model.rxns);
            model.ub(pos)  = data(j);
            match          = true;
        end
    end
end

%Count mass of non-measured enzymes:
measured       = ~isnan(model.concs);
concs_measured = model.concs(measured);
Pmeasured      = sum(concs_measured);

if Pmeasured > 0
    %Calculate fraction of non measured proteins in model out of remaining mass:    
    [fn,~] = measureAbundance(model.enzymes(~measured),'prot_abundance.txt');
    fm     = Pmeasured/Ptot;
    f      = fn/(1-fm);
    %Discount measured mass from global constrain:
    fs = (Ptot - Pmeasured)/Pbase*f*sigma;
else
    fs = f*sigma;
end

%Constrain the rest of enzymes with the pool assumption:
if sum(strcmp(model.rxns,'prot_pool_exchange')) == 0
    model = constrainPool(model,~measured,1000);
end

model = fixedModifications(model);
model = changeProtein(model,Ptot,fs);

%Display some metrics:
disp(['Total protein amount measured = '     num2str(Pmeasured)              ' g/gDW'])
disp(['Total enzymes measured = '            num2str(sum(measured))          ' enzymes'])
disp(['Enzymes in model with 0 g/gDW = '     num2str(sum(concs_measured==0)) ' enzymes'])
disp(['Total protein amount not measured = ' num2str(Ptot - Pmeasured)       ' g/gDW'])
disp(['Total enzymes not measured = '        num2str(sum(~measured))         ' enzymes'])
disp(['Total protein in model = '            num2str(Ptot)                   ' g/gDW'])

[model,enzUsages,modifications] = flexibilizeProteins(model,0.1,2);

%Plot histogram (if there are measurements):
if sum(concs_measured) > 0
    figure
    hist(concs_measured*1e3,10.^(-3:0.5:3))
    set(gca,'xscale','log')
    xlim([1e-3,1e3])
    xlabel('Protein amount [mg/gDW]');
    ylabel('Frequency');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = cleanDataset(data)
for i=1:length(data)
    if data(i)<=0
        data(i) = NaN;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%