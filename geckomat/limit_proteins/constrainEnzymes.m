 function [model,enzUsages,modifications,GAM,massCoverage] = constrainEnzymes(model,f,GAM,Ptot,pIDs,data,gRate,c_UptakeExp)
% constrainEnzymes
%
%   Main function for overlaying proteomics data on an enzyme-constrained
%   model. If chosen, also scales the protein content, optimizes GAM, and
%   flexibilizes the proteomics data.
%
%   model           ecModel.
% 	f				(Opt) Estimated mass fraction of enzymes in model.
%	GAM				(Opt) Growth-associated maintenance value. If not
%					provided, it will be fitted to chemostat data.
%   Ptot            (Opt) Total protein content, provide if desired content
%                   is different from the one reported in getModelParameters [gProt/gDw]
% 	pIDs			(Opt) Protein IDs from proteomics data.
%	data			(Opt) Protein abundances from proteomics data [mmol/gDW].
%   gRate           (Opt) Experimental growth rate at which the proteomics
%                  data were obtained [1/h]
%   c_UptakeExp     (Opt) Experimentally measured glucose uptake rate
%                   [mmol/gDW h].
%
%   model           ecModel with calibrated enzyme usage upper bounds
%   enzUsages       Calculated enzyme usages after final calibration
%                   (enzyme_i demand/enzyme_i upper bound)
%   modifications   Table with all the modified values
%                   (Protein ID/old value/Flexibilized value)
%   GAM             Fitted GAM value for the ecModel
%   massCoverage    Ratio between measured and total mass of protein in the model
%
%   Usage: [model,enzUsages,modifications, GAM,massCoverage] = constrainEnzymes(model,f,GAM,Ptot,pIDs,data,gRate,c_UptakeExp)
%
%   Benjamin J. Sanchez. Last update 2018-12-11
%   Ivan Domenzain.      Last update 2019-12-13
%

%get model parameters
cd ..
parameters = getModelParameters;
sigma      = parameters.sigma;
c_source   = parameters.c_source;
cd limit_proteins
%Compute f if not provided:
if nargin < 2
    [f,~] = measureAbundance(model.enzymes);
end
%Leave GAM empty if not provided (will be fitted later):
if nargin < 3
    GAM = [];
end
%Load Ptot if not provided:
if nargin < 4
    Ptot = parameters.Ptot;
end
%No UB will be changed if no data is available -> pool = all enzymes(FBAwMC)
if nargin < 5
    pIDs = cell(0,1);
    data = zeros(0,1);
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
    [fn,~] = measureAbundance(model.enzymes(~measured));
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
if sum(data)==0
    %Modify protein/carb content and GAM:
    [model,GAM] = scaleBioMass(model,Ptot,GAM);
end
%Display some metrics:
disp(['Total protein amount measured = '     num2str(Pmeasured)              ' g/gDW'])
disp(['Total enzymes measured = '            num2str(sum(measured))          ' enzymes'])
disp(['Enzymes in model with 0 g/gDW = '     num2str(sum(concs_measured==0)) ' enzymes'])
disp(['Total protein amount not measured = ' num2str(Ptot - Pmeasured)       ' g/gDW'])
disp(['Total enzymes not measured = '        num2str(sum(~measured))         ' enzymes'])
disp(['Total protein in model = '            num2str(Ptot)                   ' g/gDW'])
massCoverage = Pmeasured/Ptot;
enzUsages    = [];
if nargin > 7
    [model,enzUsages,modifications] = flexibilizeProteins(model,gRate,c_UptakeExp,c_source);
end

if isempty(enzUsages)
    enzUsages      = table({},zeros(0,1),'VariableNames',{'prot_IDs' 'usage'});
    modifications  = table({},zeros(0,1),zeros(0,1),'VariableNames',{'protein_IDs' 'previous_values' 'modified_values' 'flex_mass'});
else
     plotHistogram(enzUsages.usage,'Enzyme usage [-]',[0,1],'Enzyme usages','usages')
end
%Plot histogram (if there are measurements):
plotHistogram(concs_measured,'Protein amount [mg/gDW]',[1e-3,1e3],'Modelled Protein abundances','abundances')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotHistogram(variable,xlabelStr,xlimits,titleStr,option)
if iscell(variable)
    cell2mat(variable);
end
if sum(variable) > 0
    variable(variable==0) = 1E-15;
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
