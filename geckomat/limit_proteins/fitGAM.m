function GAM = fitGAM(model,verbose)
% fitGAM
%
% GAM = fitGAM(model,verbose)
% Returns a fitted GAM for the yeast model.
%
% Usage: GAM = fitGAM(model,verbose)
%
% Benjamin Sanchez. Last update: 2018-10-27
% Ivan Domenzain.   Last update: 2020-04-06

cd ..
parameters = getModelParameters;
cd limit_proteins

if isfield(parameters,'GAM')
    GAM = parameters.GAM;
else
    
    if nargin<2
        verbose = true;
    end
    if isfield(parameters,'maxGAM')
        maxGAM = parameters.maxGAM;
    else
        maxGAM = 100;
    end

    %Load chemostat data:
    fid = fopen('../../databases/chemostatData.tsv','r');
    cSource{1} = fgetl(fid);
    cSource{2} = regexprep(cSource{1},'.*\t(.*)Uptake.*','$1');
    if strcmp(cSource{1},cSource{2})
        error("Cannot recognize headerline of 'chemostatData.tsv'.");
    end
    exp_data = textscan(fid,'%f32 %f32 %f32 %f32','Delimiter','\t');
    exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4}];
    fclose(fid);
    
    %Remove limitation on enzymes (if any):
    model = setParam(model,'ub','prot_pool_exchange',+1000);
    
    %GAMs to span:
    fprintf('Estimating GAM...')
    GAM = 20:5:maxGAM;
    
    %1st iteration:
    GAM = iteration(model,GAM,exp_data,verbose);
    
    %2nd iteration:
    GAM = iteration(model,GAM-10:1:GAM+10,exp_data,verbose);
    
    %3rd iteration:
    [GAM,fitError] = iteration(model,GAM-1:0.1:GAM+1,exp_data,verbose);
    
    %If verbose output is not required, then the only displayed value is the optimal one
    if ~verbose
        fprintf(' Done!\n')
        disp(['Fitted GAM = ' num2str(GAM) ' -> Error = ' num2str(fitError)])
    end
    %Plot fit:
    mod_data = simulateChemostat(model,exp_data,GAM);
    figure
    hold on
    cols = [0,1,0;0,0,1;1,0,0];
    b    = zeros(1,length(exp_data(1,:))-1);
    for i = 1:length(exp_data(1,:))-1
        b(i) = plot(mod_data(:,1),mod_data(:,i+1),'Color',cols(i,:),'LineWidth',2);
        plot(exp_data(:,1),exp_data(:,i+1),'o','Color',cols(i,:),'MarkerFaceColor',cols(i,:))
    end
    title(['GAM fitting for growth on ', cSource{2}])
    xlabel('Dilution rate [1/h]')
    ylabel('Exchange fluxes [mmol/gDWh]')
    legend(b,[cSource{2}, ' consumption'],'O2 consumption','CO2 production','Location','northwest')
    hold off
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GAM,fitError] = iteration(model,GAM,exp_data,verbose)

fitting = ones(size(GAM))*1000;

for i = 1:length(GAM)
    %Simulate model and calculate fitting:
    mod_data   = simulateChemostat(model,exp_data,GAM(i));
    R          = (mod_data - exp_data)./exp_data;
    fitting(i) = sqrt(sum(sum(R.^2)));
    if verbose
        fprintf(['\nGAM = ' num2str(GAM(i)) ' -> Error = ' num2str(fitting(i)) '\n'])
    else
        fprintf('.')
    end
end

%Choose best:
[fitError,best] = min(fitting);

if best == 1 || best == length(GAM)
    error('GAM found is sub-optimal: please expand GAM search bounds.')
else
    GAM   = GAM(best);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mod_data = simulateChemostat(model,exp_data,GAM)
%Modify GAM withouth changing the protein content:
Pbase = sumProtein(model);
model = scaleBioMass(model,Pbase,GAM,false);
cd ..
parameters = getModelParameters;
cd limit_proteins
%Relevant positions:
exch_names  = parameters.exch_names;
pos(1)      = find(strcmp(model.rxnNames,exch_names{1}));
pos(2)      = find(strcmp(model.rxnNames,exch_names{2}));
pos(3)      = find(strcmp(model.rxnNames,exch_names{3}));
pos(4)      = find(strcmp(model.rxnNames,exch_names{4}));

%Simulate chemostats:
mod_data = zeros(size(exp_data));
for i = 1:length(exp_data(:,1))
    %Fix biomass and minimize glucose
    model = setParam(model,'lb',model.rxns(pos(1)),exp_data(i,1));
    model = setParam(model,'ub',model.rxns(pos(2)),+10);
    model = setParam(model,'obj',model.rxns(pos(2)),-1);
    sol   = solveLP(model,1);
    %Store relevant variables:
    mod_data(i,:) = sol.x(pos)';
end

end