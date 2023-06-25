function [outV,gRate] = plotCrabtree(ecModel)

gRate = [0:0.025:0.4];
outV  = zeros(numel(ecModel.rxns),numel(gRate));
glcEx = getIndexes(ecModel,'r_1714','rxns');
ecModel = setParam(ecModel,'obj','r_1714',1);
totP  = -ecModel.lb(strcmp(ecModel.rxns,'prot_pool_exchange'));

for i=1:numel(gRate)
    tmpModel = setParam(ecModel,'lb','r_2111',gRate(i));
    sol=solveLP(tmpModel);
    if ~isempty(sol.x)
        tmpModel = setParam(tmpModel,'lb','r_1714',sol.x(glcEx)*1.01);
        tmpModel = setParam(tmpModel,'obj','prot_pool_exchange',1);
        sol=solveLP(tmpModel);
        outV(:,i) = sol.x;
    end
end

% Gather experimental data
fID = fopen(fullfile(findGECKOroot,'tutorials','full_ecModel','data','vanHoek1998.tsv'),'r');
expData = textscan(fID,'%f %f %f %f %f %f %f %f','Delimiter',';','HeaderLines',2);
fclose(fID);

fluxToPlot = [expData{2} expData{3} expData{4} expData{5}];

% Get reaction indices in the model
rxnsToPlot = getIndexes(ecModel,{'r_1992','r_1672','r_1714','r_1761'},'rxns');

% Make two plots
tiledlayout(1,2)
dataColor = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
% Plot fluxes
nexttile
plot(gRate,abs(outV(rxnsToPlot,:)));
hold on
scatter(expData{1},fluxToPlot,[],dataColor)
legend(ecModel.rxnNames(rxnsToPlot),'Location','northwest')
xlabel('Growth rate (/hour)')
ylim([0 20])
ylabel('Absolute flux (mmol/gDCWh)')
hold off

% Plot total protein usage
nexttile
poolRxn = getIndexes(ecModel,'prot_pool_exchange','rxns');
plot(gRate,abs(outV(poolRxn,:))/totP)
xlabel('Growth rate (/hour)')
ylabel('Fraction of protein pool used')
fig=gcf;
fig.Position(3:4) = [700 300];
end
