function [outV,gRate] = plotCrabtree(ecModel)

gRate = [0:0.025:0.4];
outV  = zeros(numel(ecModel.rxns),numel(gRate));
glcEx = getIndexes(ecModel,'r_1714','rxns');
ecModel = setParam(ecModel,'obj','r_1714',1);

for i=1:numel(gRate)
    tmpModel = setParam(ecModel,'lb','r_2111',gRate(i));
    sol=solveLP(tmpModel);
    if ~isempty(sol.x)
        tmpModel = setParam(tmpModel,'lb','r_1714',sol.x(glcEx)*1.01);
        tmpModel = setParam(tmpModel,'obj','prot_pool_exchange',1);
        sol=solveLP(tmpModel,1);
        outV(:,i) = sol.x;
    end
end

% Gather experimental data
expData = readtable(fullfile(findGECKOroot,'tutorials','full_ecModel','data','vanHoek1998.csv'));
fluxToPlot = table2array(expData(:,2:end));

% Get reaction indices in the model
rxnsToPlot = getIndexes(ecModel,{expData.Properties.VariableNames{2:end}},'rxns');

% Make two plots
tiledlayout(1,2)
% Plot fluxes
nexttile
plot(gRate,abs(outV(rxnsToPlot,:)))
hold on
scatter(expData.r_2111,fluxToPlot)
legend(ecModel.rxnNames(rxnsToPlot))
xlabel('Growth rate (/hour)')
ylabel('Absolute flux (mmol/gDCW/h)')
hold off

% Plot total protein usage
nexttile
poolRxn = getIndexes(ecModel,'prot_pool_exchange','rxns');
plot(gRate,abs(outV(poolRxn,:))/125)
xlabel('Growth rate (/hour)')
ylabel('Fraction of protein pool used')
end
