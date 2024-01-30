% This file accompanies the GECKO 3 Nature Protocols paper https://doi.org/10.1038/s41596-023-00931-7.
%
% The function of this script is to demonstrate the reconstruction and
% analysis of a *full* ecModel. As example, it here uses the yeast-GEM
% model of Saccharomyces cerevisiae as starting point. However, this script
% does not claim to construct a "production-ready" ecYeastGEM model:
% dependent on how you intend to use the ecModel it may require additional
% curation and evaluation.
%
% DO NOT DIRECTLY USE THE ECMODEL GENERATED HERE OUTSIDE OF THIS TUTORIAL.
%
% This script has limited comments and explanations of how the various
% functions work and only the STEPs essential for this script are shown. For 
% more detailed descriptions and examples of the other STEPs, see
% tutorials/full_ecModel/protocol.m.

% Prepare software and model adapter
checkInstallation;
GECKOInstaller.install
setRavenSolver('gurobi')

% STEP 8 Load the appropriate model adapter
adapterLocation = fullfile(findGECKOroot,'tutorials','light_ecModel','HumanGEMAdapter.m');
adapter = ModelAdapterManager.setDefault(adapterLocation); 

% STEP 9 Load the starting human-GEM model
model = loadConventionalGEM();

% STEP 10-11 Reconstruct light ecModel
% Because Human-GEM contains gene identifiers in ENSEMBL format without
% trailing version identifier, and there is no field in Uniprot that
% exactly matches that format, a custom uniprotConversion.tsv is generated
% based on the genes.tsv that is distributed with Human-GEM:
% https://github.com/SysBioChalmers/Human-GEM/blob/194ebe5431c83e25f78df61caacad2fa485b5cb4/model/genes.tsv
% of which the first and the fourth column are kept to make
% data/uniprotConversion.tsv. The second parameter (true) here signifies
% that a light ecModel will be generated.
[ecModel, noUniprot] = makeEcModel(model,true);

% The above command generates a warning regarding potentially problematic
% gene associations in various grRules. This is an expected output here and
% can safely be ignored.

% noUniprot contains 6 genes that could not be mapped, so the genes.tsv
% from Human-GEM was not sufficient and additional curation should be done.
% This will be skipped for now, as 6 genes is not too many.

% STEP 12-13 Apply complex data
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

% STEP 16-19 Search for kcat values in BRENDA
ecModel         = getECfromGEM(ecModel);
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% STEP 20-25 Predict kcat values with DLKcat
ecModel         = findMetSmiles(ecModel);
% DLKcat.tsv files are ecModel-type specific and cannot be interchanged
% between full and light ecModels. This is due to a difference in how
% reactions are represented in the ecModel.ec structure (isozymatic
% reactions are split in full ecModels, and a suffix (e.g. _EXP_1) is
% attached to reaction identifiers; while in light models such reactions
% are not split and a counting prefix (e.g. 001_) is attached to reaction
% identifiers.

% The light_ecModel tutorial already comes with a DLKcat.tsv file populated
% with kcat values. If this file should be regenerated, the lines below
% should be uncommented. Note that this overwrites the existing files,
% thereby discarding existing kcat predictions.
%writeDLKcatInput(ecModel,[],[],[],[],true);
%runDLKcat();

kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 26 Merge BRENDA and DLKcat derived values
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

% STEP 27 Implement kcat values into model.ec.kcat
ecModel = selectKcatValue(ecModel,kcatList_merged);

% STEP 28
% applyCustomKcats can be run on light ecModels in a similar way as for
% full ecModels. However, no custom kcat values are defined for the ecModel
% in this tutorial.

% STEP 29
% getKcatAcrossIsozymes cannot be run on light ecModels.

% STEP 30
% getStandardKcat can be run on light ecModels in a similar way as for full
% ecModels. But due to the different model structure, no "standard"
% pseudo-enzyme is included.
% 
% Instead of introducing a "standard" pseudo-enzyme, the new
% standard enzyme usage constraints follow the light ecModel structure,
% where the prot_pool is included as pseudo-substrate with
% -standardMW/standardKcat as the stoiciometric coefficient.
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% STEP 31 Apply kcat values to S-matrix
% In light ecModels the kcat values are also kept in ecModels.ec.kcat, and
% only after running applyKcatConstraints are these introduced in the
% ecModel.S matrix. In constrast to full ecModels, the applyKcatConstraints
% function checks across isozymes which is the most efficient enzyme
% (lowest MW/kcat value), and uses this value.
ecModel = applyKcatConstraints(ecModel);

% STEP 32 Constrain with total protein pool
ecModel = setProtPoolSize(ecModel);
sol=solveLP(ecModel,1)
printFluxes(ecModel,sol.x)
saveEcModel(ecModel);

% STEP 43-44
% sensitivityTuning also works on light ecModels. However, in this tutorial
% an ecModel for the generic Human-GEM model is generated. This model is
% not directly suitable for simulations, as it does not represent a
% specific cell-type. It is therefore not trivial to set what exchange
% reactions should be open (to represent the culture medium [for cultured
% cells, or otherwise the natural environment]), and what maximum growth
% rate should be reached.

% STEP 53-57, 64-65
% Proteomics integration is NOT possible with light ecModels.

% STEP 69-70
% Simulating fluxes in light ecModels would follow a similar strategy as
% would be suitable for simulating full ecModels without proteomics
% integration. E.g. first optimization of a "classical" objective function
% (e.g. maximize growth rate, maximize ATP turnover, or minimize glucose
% uptake), constraining the ecModel with the objective value that was
% reached, followed by a second optimization where the total enzyme usage
% is minimized, yielding the most efficient enzyme allocation.
ecModel = loadEcModel();
ecModel = setParam(ecModel,'obj',adapter.params.bioRxn,1);
sol = solveLP(ecModel)
fprintf('Growth rate reached: %g /hour.\n', abs(sol.f))
% Set growth lower bound to 99% of the previous value.
ecModel = setParam(ecModel,'lb',adapter.params.bioRxn,0.99*abs(sol.f));
% Minimize protein pool usage.
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
sol = solveLP(ecModel)
fprintf('Minimum protein pool usage: %g mg/gDCW.\n', abs(sol.f))

% STEP 71
% Individual enzyme usages cannot be investigated in light ecModels, as
% these are not explicitly included in the S-matrix.

% STEP 72
% Mapping fluxes to the conventional GEM works identical as for full
% ecModels. Due to enzymes not being explicitly included in the S-matrix,
% the enzUsageFlux only includes the protein pool exchange.
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, model, sol.x);

% STEP 76-77
% The comparison between full and light ecModel is demonstrated with a
% custom function in tutorials/full_ecModel/code/plotlightVSfull.m.

% STEP 78 Make context-specific ecModel by subsetting from a general model
% To exemplify the construction of a context-specific ecModel, a
% conventional GEM of HT-29 cell line is loaded.
HT29 = readYAMLmodel(fullfile(adapter.params.path,'models','HT29-GEM.yml'));

% Make a context-specific ecModel based on the generic Human-GEM.
ecModel = loadEcModel();
ecHT29 = getSubsetEcModel(ecModel,HT29);

% STEP 79 Compare results from generic and context-specific models
% Without detailed inspection of all differences between the three models, 
% it is clear that both the enzyme constraints and contextualization have
% had its impact as the maximum growth rate is affected:
sol = solveLP(HT29);
fprintf('Growth rate in HT29-GEM: %.3f /hour.\n', abs(sol.f))
sol = solveLP(ecModel);
fprintf('Growth rate in ecHuman-GEM: %.3f /hour.\n', abs(sol.f))
sol = solveLP(ecHT29);
fprintf('Growth rate in ecHT29-GEM: %.3f /hour.\n', abs(sol.f))
