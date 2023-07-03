function [solL, solF] = plotlightVSfull()

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
adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel','YeastGEMAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);
model = loadConventionalGEM();
model = setParam(model,'lb','r_1714',-1000);
model = setParam(model,'obj','r_2111',1);

%
fprintf('Comparison of duration light vs. full ecModel\n')
tic;
ecFull   = makeEcModel(model,false);
ecFull   = getECfromGEM(ecFull);
kcatList = fuzzyKcatMatching(ecFull);
ecFull   = selectKcatValue(ecFull,kcatList);
ecFull   = applyKcatConstraints(ecFull);
ecFull   = setProtPoolSize(ecFull);
fullTime = toc;

%
tic;
ecLight  = makeEcModel(model,true);
ecLight  = getECfromGEM(ecLight);
kcatList = fuzzyKcatMatching(ecLight);
ecLight  = selectKcatValue(ecLight,kcatList);
ecLight  = applyKcatConstraints(ecLight);
ecLight  = setProtPoolSize(ecLight);
lightTime = toc;
fprintf('ecModel reconstruction: %.0f%% (%.0f vs %.0f seconds)\n', (lightTime/fullTime)*100, lightTime, fullTime);

%
tic;solFull  = solveLP(ecFull,1);
fullTime = toc;
tic;solLight = solveLP(ecLight,1);
lightTime = toc;
fprintf('FBA: %.0f%% (%.2f vs %.2f seconds)\n', (lightTime/fullTime)*100, lightTime, fullTime);

%
tic; solF = mapRxnsToConv(ecFull,model,solFull.x);
fullTime = toc;
tic; solL = mapRxnsToConv(ecLight,model,solLight.x);
lightTime = toc;
fprintf('Mapping fluxes: %.0f%% (%.3f vs %.3f seconds)\n', (lightTime/fullTime)*100, lightTime, fullTime);

%
fprintf('Growth rate that is reached: %.4f vs %.4f\n', abs(solFull.f) , abs(solLight.f))

%
scatter(abs(solF),abs(solL))
xlabel('Full ecModel')
ylabel('Light ecModel')
title('Absolute fluxes (mmol/gDCWh)')
set(gca,'yscale','log','YminorTick','on')
set(gca,'xscale','log','XminorTick','on')
set(gca,'FontSize',11)
text(1e-7,3,'Growth rate(/hour)','FontSize',11)
text(1e-7,0.5,'full:','FontSize',11)
text(1.e-6,0.5,num2str(abs(solFull.f)),'FontSize',11)
text(1e-7,0.1,'light:','FontSize',11)
text(1.e-6,0.1,num2str(abs(solLight.f)),'FontSize',11)
saveas(gca, fullfile(findGECKOroot,'tutorials','full_ecModel','output','lightVSfull.pdf'))
