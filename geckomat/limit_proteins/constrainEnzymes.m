%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = constrainEnzymes(model,Ptot,sigma,f,GAM,pIDs,data,gRate,GlucUptake)
% 
% Benjamin J. Sanchez. Last edited: 2018-10-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,enzUsages,modifications] = constrainEnzymes(model,Ptot,sigma,f,GAM,pIDs,data,gRate,GlucUptake)

%Compute f if not provided:
if nargin < 4
    [f,~] = measureAbundance(ecModel.enzymes);
end

%Leave GAM empty if not provided (will be fitted later):
if nargin < 5
    GAM = [];
end

%No UB will be changed if no data is available -> pool = all enzymes(FBAwMC)
if nargin < 6
    pIDs          = cell(0,1);
    data          = zeros(0,1);
    enzUsages     = zeros(0,1);
    modifications = cell(0,1);
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

%Get protein content in biomass pseudoreaction:
Pbase = sumProtein(model);

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
    model = constrainPool(model,~measured,full(fs*Pbase));
end

%Modify protein/carb content and GAM:
model = scaleBioMass(model,Ptot,GAM);

%Display some metrics:
disp(['Total protein amount measured = '     num2str(Pmeasured)              ' g/gDW'])
disp(['Total enzymes measured = '            num2str(sum(measured))          ' enzymes'])
disp(['Enzymes in model with 0 g/gDW = '     num2str(sum(concs_measured==0)) ' enzymes'])
disp(['Total protein amount not measured = ' num2str(Ptot - Pmeasured)       ' g/gDW'])
disp(['Total enzymes not measured = '        num2str(sum(~measured))         ' enzymes'])
disp(['Total protein in model = '            num2str(Ptot)                   ' g/gDW'])

if nargin > 7
    [model,enzUsages,modifications] = flexibilizeProteins(model,gRate,GlucUptake);
    plotHistogram(enzUsages,'Enzyme usage [-]',[0,1],'Enzyme usages','usages')
end

%Plot histogram (if there are measurements):
plotHistogram(concs_measured,'Protein amount [mg/gDW]',[1e-3,1e3],'Modelled Protein abundances','abundances')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotHistogram(variable,xlabelStr,xlimits,titleStr,option)
if sum(variable) > 0
    figure
    if strcmpi(option,'abundances')
        hist(variable*1e3,10.^(-3:0.5:3))
        set(gca,'xscale','log')
    else
        hist(variable,(0:0.05:1))
    end
    xlim(xlimits)
    xlabel(xlabelStr)
    ylabel('Frequency');
    title(titleStr)
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
