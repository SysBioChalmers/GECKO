% For reference, code is provided here that was used to make a cell-line
% specific (HT-29) model of human-GEM with the ftINIT function of RAVEN.
% This model, stored at tutorials/light_ecModel/models/HT-29-GEM.yml, is
% for demontration purposes only. The intention of the code in this file is
% not to routinely regenerate this model file, but merely as reference how
% the HT-29-GEM was constructed.

% Have the Human-GEM repo at release 1.15.0, and added to the MATLAB path.
humanGEMroot = HumanGEMInstaller.getHumanGEMMainPath;
load(fullfile(humanGEMroot,'model','Human-GEM.mat');
m = ihuman;
m.grRules = simplifyGrRules(m.grRules, true);

prepDataHumanGEM = prepHumanModelForftINIT(m, true, ...
                    fullfile(humanGEMroot,'data','metabolicTasks','metabolicTasks_Essential.txt'), ...
                    fullfile(humanGEMroot,'model','reactions.tsv'));
save('prepDataHumanGEMGeneSymbols.mat','prepDataHumanGEM');

% Use RNA-seq data of cell line HT-29 (ACH-000552) from DepMap, downloaded
% from here: https://depmap.org/portal/download/all/.

tbl = readtable('OmicsExpressionProteinCodingGenesTPMLogp1.csv');
% Cell lines are per row, genes are columns.
sel = strcmp(tbl.Var1, 'ACH-000552');
genes = tbl.Properties.VariableNames(2:end);
for x = 1:length(genes)
     a = split(genes{x}, '_');
     genes{x} = a{1};
end

% Invert Logp1
dataTmp = table2array(tbl(sel,2:end));
data = 2.^dataTmp - 1;

arrayData = struct(); 
arrayData.genes = genes.';
arrayData.levels = data.';
arrayData.tissues = {'1'};
arrayData.threshold = 1;
HT29 = ftINIT(prepDataHumanGEM,arrayData.tissues{1},[],[],arrayData,{},getHumanGEMINITSteps('1+1'),false,true,[]);

% Restore ENSEMBL gene identifiers and remove unnecessary fields
[HT29.grRules,HT29.genes,HT29.rxnGeneMat] = translateGrRules(HT29.grRules,'ENSG','Symbol');
HT29 = rmfield(HT29,{'geneShortNames','rxnReferences','rxnFrom','metFrom','rxnConfidenceScores','rxnNotes','metCharges','inchis','metFormulas','subSystems','eccodes'});
HT29.name = 'GEM of HT-29 cell-line, for use in GECKO3 tutorial';
writeYAMLmodel(HT29,fullfile('..','models','HT29-GEM.yml'));
