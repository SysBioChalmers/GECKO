function resultStruct = diff_FluxDist_analysis(refModel,condModel,bioRxn,Csource,Drate,fileName,pathwayStr)
%diff_FluxDist_analysis
%
% Function that runs parsimonious FBA simulations (specific for ecModels)
% on two condition-specific ecModels. Differential analysis is performed
% at the levels of individual metabolic reactions, enzymes (also enzyme
% usages, if proteomics data were incorporated) and subSystems/pathways.
%
% The models used for this analysis should already come with their
% condition-specific constraints.
%
%   refModel    An ecModel structure for the reference condition
%   condModel   An ecModel structure for the specific condition to test
%   bioRxn      Rxn ID for the biomass production reaction
%   Csource     Rxn name for the main carbon source uptake reaction
%               'D-glucose exchange (reversible)' in the case of glucose.
%   Drate       Dilution rate used in the chemostat experiments [1/h]. If
%               batch conditions are desired for fluxDist comparison then 
%               provide an empty vector.
%   fileName    String for the root of the file name for the results output
%               Example: 'Kma_diffFluxDist_ref_HiT'
%   pathwayStr  String indicating a pathway name if specific output files
%               for this are also desired.
%
%   resultStruct
%       rxns        Table for the differential flux analysis at the reaction
%                   level including columns for: rxnIds/ref_FluxValue/
%                   condition_FluxValue/Log2FC/metabolis subsystems for the reaction/
%                   genes associated with each reaction. 
%       proteins    Table for the differential flux analysis at the enzymes
%                   level including columns for: uniprotId/ref_usageValue (mmol/gDw h)/
%                   cond_usageValue (mmol/gDw h)/Log2FC/usagesRef/UsagesCond 
%                   (percentage usages if proteomics measurement were incorporated 
%                   for that specific protein)/metabolis subsystems for the protein/
%                   genes associated with each protein.  
%
%       If the pathwayStr is provided, then subset tables are also created 
%       at the rxn and protein level for the specified pathway.
%
%       All created tables are also stored as tab separated .txt files in
%       the 'Results' subfolder of this repository.
%
%   Usage: resultStruct = diff_FluxDist_analysis(refModel,condModel,bioRxn,Csource,Drate,fileName,pathwayStr)
%
% Last modified.  Ivan Domenzain 2019-06-10

pos(1)    = find(strcmpi(refModel.rxnNames,Csource));
pos(2)    = find(strcmpi(refModel.rxns,bioRxn));
%Get pFBA solutions for both models
if ~isempty(Drate)
    ref_sol   = simulateChemostat(refModel,Drate,pos,true);
    cond_sol  = simulateChemostat(condModel,Drate,pos,true);
else
    ref_sol   = ecModel_batchSimulation(refModel,pos(2));
    cond_sol  = ecModel_batchSimulation(condModel,pos(2));
end
%Get metabolic rxns indexes for both models
metIndxs_ref   = find(~contains(refModel.rxnNames,'prot_'));
metIndxs_cond  = find(~contains(condModel.rxnNames,'prot_'));
if isempty(setdiff(metIndxs_ref,metIndxs_cond))
    %Get enzymes usages values
    protsIndxs_ref  = [];
    protsIndxs_cond = [];
    %Proteins are not necessarily measured in both conditions, therefore
    %usage pseudoreactions may not be in the same positions when comparing
    %models
    for i=1:length(refModel.enzymes)
        enzyme          = refModel.enzymes{i};
        index           = find(contains(refModel.rxnNames,['prot_' enzyme]));
        protsIndxs_ref  = [protsIndxs_ref; index];
        index           = find(contains(condModel.rxnNames,['prot_' enzyme]));
        protsIndxs_cond = [protsIndxs_cond; index];
    end
else
    disp('Inconsistent mapping')
end
mkdir('../Results/diff_FluxDist_analyisis')
fileName = ['../Results/diff_FluxDist_analyisis/' fileName];
%Extract values from flux distributions
metFluxes_ref    = ref_sol(metIndxs_ref);
metFluxes_cond   = cond_sol(metIndxs_cond);
comp_rxns_Table  = comparativeTable(refModel,condModel,metIndxs_ref,metIndxs_cond,metFluxes_ref,metFluxes_cond);
fileNameStr      = [fileName '_rxns.txt'];
writetable(comp_rxns_Table,fileNameStr,'Delimiter','\t','QuoteStrings',false)

%Extract values from prot usage distributions
protFlux_ref    = ref_sol(protsIndxs_ref);
protFlux_cond   = cond_sol(protsIndxs_cond);
comp_prot_Table = comparativeTable(refModel,condModel,protsIndxs_ref,protsIndxs_cond,protFlux_ref,protFlux_cond,true);
fileNameStr     = [fileName '_enzymes.txt'];
writetable(comp_prot_Table,fileNameStr,'Delimiter','\t','QuoteStrings',false)

%Output structure
resultStruct.proteins = comp_prot_Table;
resultStruct.rxns     = comp_rxns_Table;
%If a subsystem/pathway of interest was provided then a subset of the tabl
%econtaining just the changes for it will then be created.
if nargin>6
    subsetTableProt = extractSubSystemFromTable(comp_prot_Table,pathwayStr);
    subsetTableRxns = extractSubSystemFromTable(comp_rxns_Table,pathwayStr);
    resultStruct.Pathway_rxns = subsetTableRxns;
    resultStruct.Pathway_prot = subsetTableProt;
    fileNameStr      = [fileName '_' pathwayStr '_rxns.txt'];
    writetable(subsetTableRxns,fileNameStr,'Delimiter','\t','QuoteStrings',false)
    fileNameStr      = [fileName '_' pathwayStr '_enzymes.txt'];
    writetable(subsetTableProt,fileNameStr,'Delimiter','\t','QuoteStrings',false)
end

end
%--------------------------------------------------------------------------
function compTable = comparativeTable(refModel,condModel,refIndxs,condIndxs,refVals,condVals,protFlag)
if nargin<7 
    protFlag = false;
end
%Extract values from flux distributions
rxnIDs_ref  = refModel.rxnNames(refIndxs);
rxnIDs_cond = condModel.rxnNames(condIndxs);
if protFlag 
    %Extract protein IDs from usage reactions
    rxnIDs_ref  = getProteinID(rxnIDs_ref);
    rxnIDs_cond = getProteinID(rxnIDs_cond);
end
subSystems     = {};
log2FC         = {};
nonZero_ref    = refVals;
nonZero_ref(nonZero_ref==0) = 1E-12;
nonZero_cond  = condVals;
nonZero_cond(nonZero_cond==0) = 1E-12;
if isempty(setdiff(rxnIDs_ref,rxnIDs_cond))
    for i=1:length(refIndxs)
        index = refIndxs(i);
        %Get subsystems for each reaction
        if isempty(refModel.subSystems{index})
            str = ' ';
        else
        	str = strjoin(refModel.subSystems{index},' // ');
        end
        subSystems = [subSystems; {str}];
        %Get fluxes FC
        Fchange = log2(nonZero_cond(i)/nonZero_ref(i));
        log2FC  = [log2FC; Fchange];
    end
    grRules        = refModel.grRules(refIndxs);
    modelIDs       = refModel.rxns(refIndxs);
    if protFlag
        %Calculate percentage enzyme usages for measured enzymes and
        %append it to the results comparative table
        usagesRef  = getEnzUsagesPercentage(refModel,rxnIDs_ref,refVals);
        usagesCond = getEnzUsagesPercentage(condModel,rxnIDs_cond,condVals);
        subSystems = mapEnzymeSubSystems(rxnIDs_cond,refModel);
        compTable  = table(rxnIDs_ref,refVals,condVals,log2FC,usagesRef,usagesCond,subSystems,grRules);
    else
        compTable  = table(rxnIDs_ref,modelIDs,refVals,condVals,log2FC,subSystems,grRules);
    end
end
end
%--------------------------------------------------------------------------
function prots = getProteinID(rxnIDs)
    rxnIDs = strrep(rxnIDs,'prot_','');
    rxnIDs = strrep(rxnIDs,'draw_','');
    prots  = strrep(rxnIDs,'_exchange','');
end
%--------------------------------------------------------------------------
function usages = getEnzUsagesPercentage(model,protIDs,values)
usages = NaN(length(protIDs),1);
for i=1:length(protIDs)
    protein = protIDs{i};
    rxnPos  = find(contains(model.rxns,['prot_' protein]));
    UB      = model.ub(rxnPos);
    if UB<1000
        usage     = values(i)/UB;
        usages(i) = usage;
    end
end
end
%--------------------------------------------------------------------------
function subsetTable = extractSubSystemFromTable(compTable,pathwayStr)
subSystems  = compTable.subSystems;
subSystems  = lower(subSystems);
indexes     = find(contains(subSystems,pathwayStr));
subsetTable = compTable(indexes,:);
end
%--------------------------------------------------------------------------
function solVector = ecModel_batchSimulation(model,pos)
model.c(:)   = 0;
model.c(pos) = 1;
protPos      = find(contains(model.rxnNames,'prot_') & contains(model.rxnNames,'_exchange'));
solVector   = zeros(length(model.rxns),1);
sol = solveLP(model);
if ~isempty(sol.x)
    model.lb(pos)    = 0.9999*sol.x(pos);
    model.c(:)       = 0;
    model.c(protPos) = -1;
    sol              = solveLP(model,1);
    if ~isempty(sol.x)
        solVector = sol.x;
    end
end
end










