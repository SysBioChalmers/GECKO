%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [GAM,fitting,results] = fitGAM(model)
% 
% Benjamin Sanchez. Last update: 2018-08-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GAM = fitGAM(model)

%Load chemostat data:
fid = fopen('../../databases/chemostatData.tsv','r');
exp_data = textscan(fid,'%f32 %f32 %f32 %f32','Delimiter','\t','HeaderLines',1);
exp_data = [exp_data{1} exp_data{2} exp_data{3} exp_data{4}];
fclose(fid);

%Remove limitation on enzymes (if any):
model = setParam(model,'ub','prot_pool_exchange',+1000);

%GAMs to span:
disp('Estimating GAM:')
GAM = 20:5:50;

%1st iteration:
GAM = iteration(model,GAM,exp_data);

%2nd iteration:
GAM = iteration(model,GAM-10:1:GAM+10,exp_data);

%3rd iteration:
GAM = iteration(model,GAM-1:0.1:GAM+1,exp_data);

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
title('GAM fitting for growth on glucose minimal media')
xlabel('Dilution rate [1/h]')
ylabel('Exchange fluxes [mmol/gDWh]')
legend(b,'Glucose consumption','O2 consumption','CO2 production','Location','northwest')
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GAM = iteration(model,GAM,exp_data)

fitting = ones(size(GAM))*1000;

for i = 1:length(GAM)
    %Simulate model and calculate fitting:
    mod_data   = simulateChemostat(model,exp_data,GAM(i));
    R          = (mod_data - exp_data)./exp_data;
    fitting(i) = sqrt(sum(sum(R.^2)));
    disp(['GAM = ' num2str(GAM(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[~,best] = min(fitting);

if best == 1 || best == length(GAM)
    error('GAM found is sub-optimal: please expand GAM search bounds.')
else
    GAM = GAM(best);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mod_data = simulateChemostat(model,exp_data,GAM)

%Modify GAM withouth changing the protein content:
[~,Pbase,~,~,~,~] = sumBioMass(model);
model = scaleBioMass(model,Pbase,GAM,false);

%Relevant positions:
pos(1) = find(strcmp(model.rxnNames,'growth'));
pos(2) = find(strcmp(model.rxnNames,'D-glucose exchange (reversible)'));
pos(3) = find(strcmp(model.rxnNames,'oxygen exchange (reversible)'));
pos(4) = find(strcmp(model.rxnNames,'carbon dioxide exchange'));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
