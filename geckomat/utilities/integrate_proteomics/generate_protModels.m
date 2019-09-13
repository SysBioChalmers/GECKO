function exchFlux_errors = generate_protModels(ecModel,grouping,flexFactor,oxPhosIDs,ecModel_batch,protBasis)
% generate_protModels
%
% Function that takes an ecModel and constraints it with absolute proteomics 
% datasets for the experimental conditions provided by the user.
%
%   ecModel       (structure) ecModel MATLAB structure (without the total protein pool constraint)
%   grouping      (vector) Number of biological replicates for each
%                 experimental condition in the dataset.
%   flexFactor    (double- optional) Flexibilization factor for carbon source 
%                 uptake flux for automated flexibilization of protein levels.
%   modelName     (string) Model name string (e.g. ecYeastGEM for S. cerevisiae)
%   oxPhosIDs     (cell) IDs for the oxidative phosphorylation related. The
%                 mean expression level across subunits is taken as the 
%                 abundance levelapplied to each of the enzyme complexes 
%                 in these reactions to avoid overconstraining.
%   ecModel_batch (structure) ecModel MATLAB structure with total protein pool
%                 constraint, if provided the process of fitting GAM for
%                 each new protein content is speed-up.
%   protBasis     (logical, default false) TRUE if proteomics dataset is
%                 provided in units of [mmol/g protein]. Default units are:
%                 [mmol/gDw]
%
% Usage:  exchFlux_errors = generate_protModels(ecModel,grouping,flexFactor,oxPhosIDs,ecModel_batch,protBasis)
%
% Last modified.  Ivan Domenzain 2019-09-13

close all
current = pwd;

%This funcion allows for flexibilization of protein absolute abundances in 
%case that ecModelP is not feasible using the automatically flexibilized 
%data, if flex factor is not specified then a factor of 1 is assumed.
if nargin<6
    protBasis = false;
    if nargin<5
        ecModel_batch = [];
        if nargin<4
            oxPhosIDs = [];
            if nargin<3
                flexFactor = 1.05;
            end
        end
    end
end
%get model parameters
cd ../..
parameters = getModelParameters;
Ptot_model = parameters.Ptot;
c_source   = parameters.c_source;
bioRXN     = parameters.bioRxn;
NGAM       = parameters.NGAM;
exch_ids   = parameters.exch_names(2:end);
%create subfolder for ecModelProts output files
mkdir('../models/prot_constrained')
%Get oxPhos related rxn IDs
oxPhos     = getOxPhosRxnIDs(ecModel,oxPhosIDs);

%Get indexes for carbon source uptake and biomass pseudoreactions
positionsEC(1) = find(strcmpi(ecModel.rxnNames,c_source));
positionsEC(2) = find(strcmpi(ecModel.rxns,bioRXN));
%Remove and substitute files and databases in GECKO
removeFile('../Databases/relative_proteomics.txt')
%Remove prot_abundance.txt  and relative_proteomics.txt files
%(for f factor calculation)
removeFile('../Databases/prot_abundance.txt')

%Load absolute proteomics dataset [mmol/gDw]
%and fermentation data (GUR, OUR, CO2 production, byProducts, Ptot, Drate)
cd utilities/integrate_proteomics/
[uniprotIDs,absValues,fermData,byProducts] = load_Prot_Ferm_Data(grouping);
conditions = fermData.conds;
Ptot       = fermData.Ptot;
Drate      = fermData.Drate;
GUR        = fermData.GUR;
CO2prod    = fermData.CO2prod;
OxyUptake  = fermData.OxyUptake;
byP_flux   = fermData.byP_flux;
exch_error = zeros(length(conditions),1);
%For each condition create a protein constrained model
for i=1:length(conditions)
    cd (current)
    disp(conditions{i})
    %Extract data for the i-th condition
    abundances   = cell2mat(absValues(1:grouping(i)));
    initialProts = uniprotIDs;
    if protBasis %If dataset units are [mmol/g prot] this converts them to: [mmol/gDw]
        abundances = abundances*Ptot(i);
    end
    absValues  = absValues(grouping(i)+1:end);
    %Filter data
    [pIDs, abundances] = filter_ProtData(uniprotIDs,abundances,1.96,true);
    filteredProts      = pIDs;
    cd ..
    for j=1:length(oxPhos)
        [abundances,pIDs] = fixComplex(oxPhos{j},ecModel,abundances,pIDs);
    end

    %Set minimal medium
    cd ../kcat_sensitivity_analysis
    ecModelP = changeMedia_batch(ecModel,c_source);
    
    %If the relative difference between the ecModel's protein content and
    %the Ptot for i-th condition is higher than 5% then biomass should be
    %rescaled and GAM refitted to this condition.
    %For fitting GAM a functional model is needed therefore an ecModel with
    %total protein pool constraint should be used
    Prot_diff = abs(Ptot_model-Ptot(i))/Ptot_model;
    if Prot_diff>=0.05
        if ~isempty(ecModel_batch)
            tempModel = changeMedia_batch(ecModel_batch,c_source);
            cd ../limit_proteins
            %scaleBiomass gets a GAM value fitted for the provided biomass
            %composition
            [~,GAM] = scaleBioMass(tempModel,Ptot(i),[],true);
        else
            cd ../limit_proteins
            [~,~,~,GAM] = constrainEnzymes(ecModelP,0.5,[],Ptot(i));
        end
        %Then the GAM and new biomass composition are set in ecModelP, which 
        %is not functional yet but should be used for incorporation of 
        %proteomics data
        ecModelP = scaleBioMass(ecModelP,Ptot(i),GAM,true);
        disp(' ')
    end
    %Block production of non-observed metabolites before data incorporation
    %and flexibilization
    expData  = [GUR(i),CO2prod(i),OxyUptake(i)];
    ecModelP = DataConstrains(ecModelP,byProducts,byP_flux(i,:));

    %Reescale biomass composition (according to the provided Ptot) and fit
    %GAM
    %Incorporate protein abundances into ecModel
    f       = 0.5; %Protein mass in model/Total theoretical proteome
    flexGUR = flexFactor*GUR(i);
    disp(['Incorporation of proteomics constraints for ' conditions{i} ' condition'])
    [ecModelP,usagesT,modificationsT,~,coverage] = constrainEnzymes(ecModelP,f,GAM,Ptot(i),pIDs,abundances,Drate(i),flexGUR);
    matchedProteins = usagesT.prot_IDs;
    disp(' ')
    disp(['The total number of proteins in the dataset is:                ' num2str(length(initialProts))])
    disp(['The total number of proteins in the filtered dataset is:       ' num2str(length(filteredProts))])
    disp(['The total number of filtered proteins present in the model is: ' num2str(length(matchedProteins))])
    disp(['The mass ratio between measured and unmeasured protein is:     ' num2str(coverage)])
    writeProtCounts(initialProts,filteredProts,matchedProteins,ecModelP.enzymes,coverage); 
    %Set chemostat conditions constraints and fit NGAM
    cd (current)
    %NGAM interval for fitting
    interval = [0 5];
    ecModelP = setStressConditions(ecModelP,Drate(i),positionsEC,expData,NGAM,interval);
    
    %Get optimal flux distribution and display exchange fluxes
    solution = solveLP(ecModelP,1);
    if ~isempty(solution.f)
        printFluxes(ecModelP,solution.x,true,1E-4)
    end
    %Fix experimental Glucose uptake rate and save models
    cd ..
    ecModelP      = setChemostatConstraints(ecModelP,positionsEC,Drate(i),true,0.01,GUR(i));
    exch_error(i) = getExchanges(ecModelP,exch_ids,expData);
    save(['../../models/prot_constrained/ecModel_Prot_' conditions{i} '.mat'],'ecModelP')
    %save .txt file
    writetable(usagesT,['../../models/prot_constrained/enzymeUsages_' conditions{i} '.txt'],'Delimiter','\t')
    writetable(modificationsT,['../../models/prot_constrained/modifiedEnzymes_' conditions{i} '.txt'],'Delimiter','\t')
end
exchFlux_errors = table(conditions,num2cell(exch_error),'VariableNames',{'conditions' 'avg_error'});
writetable(exchFlux_errors,'../../models/prot_constrained/exchangeFluxes_avgError.txt','Delimiter','\t')
cd (current)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function condModel = setStressConditions(model,Drate,pos,expData,NGAM,interval)
minProt = true;
ecFlag  = true;
%Rescale biomass set fitted GAM for std conditions and fit non-growth
%associated maintenance for stress conditions
cd ..
condModel = setChemostatConstraints(model,pos,Drate,minProt,0.01);
cd integrate_proteomics
disp(' ')
condModel = fitNGAM(condModel,NGAM,expData,interval,ecFlag);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function removeFile(fileName)
if exist(fileName,'file')~= 0
    delete(fileName)
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = DataConstrains(model,byProducts,bounds,fixed)
if nargin <4
    fixed = false;
end
if ~isempty(byProducts)
    disp('Constraining byproducts exchange fluxes with fermentation data')
    for i=1:length(byProducts)
        %Get exchange rxn index
        rxnName = [byProducts{i} ' exchange'];
        BPindex = find(strcmpi(model.rxnNames,rxnName));
        if ~isempty(BPindex)
            disp([byProducts{i} ' exchange has been constrained to: ' num2str(bounds(i)) ' [mmol/gDw h]'])
            %Allow a 5% of flexibility 
            model = setParam(model,'ub',BPindex,1.05*bounds(i));
            if fixed
                model = setParam(model,'lb',BPindex,0.95*bounds(i));
            end
        else
            disp(['No exchange rxn for ' byProducts{i} ' was found in ecModel'])
        end
    end
end
disp(' ')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function oxPhos = getOxPhosRxnIDs(model,rxnIDs)
if isempty(rxnIDs) %S. cerevisiae case (yeast8) 
    rxnIDs = {'r_1021' 'r_0439' 'r_0438' 'r_0226'}; 
end    
oxPhos = [];
for i=1:length(rxnIDs)
    ID = rxnIDs{i};
    isoEnzymes = model.rxns(find(contains(model.rxns,ID)));
    isoEnzymes = isoEnzymes(~contains(isoEnzymes,'arm_'));
    oxPhos = [oxPhos; isoEnzymes];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function avg_error = getExchanges(model,exch_ids,expData)
exch_ids  = [exch_ids(1) exch_ids(3) exch_ids(2)];
for i=1:length(exch_ids)
    pos(i) = find(strcmp(model.rxnNames,exch_ids(i)));
end
sol = solveLP(model,1);
if ~isempty(sol.x)
    exchanges = sol.x(pos)';
else
    exchanges = zeros(1,length(length(pos)));
end
residues  = (abs(exchanges) - expData)./expData;
avg_error = mean(abs(residues));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeProtCounts(initial_prots,filtered_prots,matched_prots,model_prots,mass_coverage)
initial_prots  = length(initial_prots);
filtered_prots = length(filtered_prots);
matched_prots  = length(matched_prots);
model_prots    = length(model_prots);
T = table(initial_prots,filtered_prots,matched_prots,model_prots,mass_coverage);
writetable(T,'../../models/prot_constrained/prot_counts.txt','Delimiter','\t')
end
