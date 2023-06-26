function plotlightVSfull()

%% Compare light with full ecModel. 
% 
% 1. To minimize the difference between the two types of ecModels beyond their
% formulation, two relatively simple ecModels (only use fuzzyKcatMatching
% as kcat source, no kcat tuning etc.) are generated for S. cerevisiae.
% 2. Maximum growth rate is simulated in both ecModels, unconstrained by
% glucose uptake.
% 3. The flux distributions are mapped back to the conventional GEM to
% allow for comparison.
% 4. A plot is made and stored in ../output/lightVSfull.pdf.

%
disp('Prepare GECKO adapter and load conventional GEM')
adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel','YeastGEMAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);
model = loadConventionalGEM();
model = setParam(model,'lb','r_1714',-1000);
model = setParam(model,'obj','r_2111',1);

%
disp('Make full ecModel')
ecFull   = makeEcModel(model,false);
ecFull   = getECfromGEM(ecFull);
kcatList = fuzzyKcatMatching(ecFull);
ecFull   = selectKcatValue(ecFull,kcatList);
ecFull   = applyKcatConstraints(ecFull);
ecFull   = setProtPoolSize(ecFull);

%
disp('Make light ecModel')
ecLight  = makeEcModel(model,true);
ecLight  = getECfromGEM(ecLight);
kcatList = fuzzyKcatMatching(ecLight);
ecLight  = selectKcatValue(ecLight,kcatList);
ecLight  = applyKcatConstraints(ecLight);
ecLight  = setProtPoolSize(ecLight);

%
disp('Perform FBA, maximizing growth')
solFull  = solveLP(ecFull,1);
solLight = solveLP(ecLight,1);

%
disp('Growth rate that is reached:')
fprintf('full ecModel:  %f \nlight ecModel: %f\n', abs(solFull.f) , abs(solLight.f))

%
disp('Map reaction rates back to conventional GEM')
solFx = mapRxnsToConv(ecFull,model,solFull.x);
solLx = mapRxnsToConv(ecLight,model,solLight.x);

%
disp('Prepare and save plot')
scatter(abs(solFx),abs(solLx))
xlabel('Full ecModel')
ylabel('Light ecModel')
title('Absolute fluxes (mmol/gDCWh)')
set(gca,'yscale','log','YminorTick','on')
set(gca,'xscale','log','XminorTick','on')
set(gca,'FontSize',11)
text(1e-9,3,'Growth rate(/hour)','FontSize',11)
text(1e-9,0.5,'full:','FontSize',11)
text(1.e-8,0.5,num2str(abs(solFull.f)),'FontSize',11)
text(1e-9,0.1,'light:','FontSize',11)
text(1.e-8,0.1,num2str(abs(solLight.f)),'FontSize',11)
saveas(gca, fullfile('..','output','lightVSfull.pdf'))
