function generate_protModels(ecModel,grouping,name,flexFactor,ecModel_batch)
% generate_protModels
%
% Function that takes an ecModel and constraints it with absolute proteomics 
% datasets for the experimental conditions provided by the user.
%
%   ecModel       (structure) ecModel MATLAB structure (without the total protein pool constraint)
%   grouping      (vector) Number of biological replicates for each
%                 experimental condition in the dataset.
%   name          (string) name string for the ecModel + proteomics and its
%                 container folder (GECKO/models/prot_constrained/name/name...)
%   flexFactor    (double- optional) Flexibilization factor for carbon source 
%                 uptake flux for automated flexibilization of protein levels.
%   ecModel_batch (structure) ecModel MATLAB structure with total protein pool
%                 constraint, if provided the process of fitting GAM for
%                 each new protein content is speed-up.
%
% Usage:  generate_protModels(ecModel,grouping,name,flexFactor,oxPhosIDs,ecModel_batch)
%
% Last modified.  Ivan Domenzain 2020-01-20
close all
current = pwd;

%This funcion allows for flexibilization of protein absolute abundances in 
%case that ecModelP is not feasible using the automatically flexibilized 
%data, if flex factor is not specified then a factor of 1 is assumed.
if nargin<5
    ecModel_batch = [];
    if nargin<4
        flexFactor = 1.05;
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
%Get oxPhos related rxn IDs
oxPhos     = getOxPhosRxnIDs(ecModel,parameters);
%create subfolder for ecModelProts output files
mkdir(['../models/prot_constrained/' name])
%Get indexes for carbon source uptake and biomass pseudoreactions
positionsEC(1) = find(strcmpi(ecModel.rxnNames,c_source));
positionsEC(2) = find(strcmpi(ecModel.rxns,bioRXN));
%Remove prot_abundance.txt  and relative_proteomics.txt files
%(for f factor calculation)
try
    movefile ../Databases/prot_abundance.txt ../Databases/prot_abundance_temp.txt
catch
    disp('prot_abundance.txt file not found in Databases folder') 
end
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
error_RQ   = zeros(length(conditions),1);
RQ         = zeros(length(conditions),1);
%For each condition create a protein constrained model
for i=1:length(conditions)
    cd (current)
    disp(conditions{i})
    %Extract data for the i-th condition
    abundances   = cell2mat(absValues(1:grouping(i)));
    initialProts = uniprotIDs;
    absValues    = absValues(grouping(i)+1:end);
    %Filter data
    [pIDs, abundances] = filter_ProtData(uniprotIDs,abundances,1.96,true);
    filteredProts      = pIDs;
    cd ..
    %correct oxPhos complexes abundances
    if ~isempty(oxPhos)
        for j=1:length(oxPhos)
            [abundances,pIDs] = fixComplex(oxPhos{j},ecModel,abundances,pIDs);
        end
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
    ecModelP = DataConstrains(ecModelP,byProducts,byP_flux(i,:),1.1);
    %Constrain CO2 production and oxygen uptake to minimize weird exchanges
    ecModelP = DataConstrains(ecModelP,{'carbon dioxide'},expData(2),[Inf 0.9]);
    
    %Reescale biomass composition (according to the provided Ptot) and fit
    %GAM
    %Incorporate protein abundances into ecModel
    f       = 0.5; %Protein mass in model/Total theoretical proteome
    flexGUR = flexFactor*GUR(i);
    disp(['Incorporation of proteomics constraints for ' conditions{i} ' condition'])
    [ecModelP,usagesT,modificationsT,~,coverage] = constrainEnzymes(ecModelP,f,GAM,Ptot(i),pIDs,abundances,Drate(i),flexGUR);
    matchedProteins = usagesT.prot_IDs;
    disp(' ')
    prot_input = {initialProts filteredProts matchedProteins ecModelP.enzymes coverage};
    writeProtCounts(conditions{i},prot_input,name); 
    %Set chemostat conditions constraints and fit NGAM
    cd (current)
    %NGAM interval for fitting
    interval = [0 5];
    ecModelP = setStressConditions(ecModelP,Drate(i),positionsEC,expData,NGAM,interval);
    %Fix experimental Glucose uptake rate and save models
    cd ..
    ecModelP = setChemostatConstraints(ecModelP,positionsEC,Drate(i),true,0.01,GUR(i));
    %Get optimal flux distribution and display exchange fluxes
    solution = solveLP(ecModelP,1);
    if ~isempty(solution.f)
        fileFluxes = ['../../models/prot_constrained/' name '/fluxes_Exch_' conditions{i} '.txt'];
        printFluxes(ecModelP,solution.x,true,1E-4,fileFluxes)
        %If model is feasible then compute RQ and compare to experimental
        %data
        [RQ(i), error_RQ(i)]= getRQ(ecModelP,solution.x,exch_ids,expData);
        disp(['The error in the respiratory quotient prediction is: ' num2str(error_RQ(i)*100) '%'])
    end
    save(['../../models/prot_constrained/' name '/' name '_' conditions{i} '.mat'],'ecModelP')
    %save .txt file
    writetable(usagesT,['../../models/prot_constrained/' name '/enzymeUsages_' conditions{i} '.txt'],'Delimiter','\t')
    writetable(modificationsT,['../../models/prot_constrained/' name '/modifiedEnzymes_' conditions{i} '.txt'],'Delimiter','\t')
end
RQ_predictions = table(conditions,num2cell(RQ),num2cell(error_RQ),'VariableNames',{'conditions' 'respiratory_quotient' 'error'});
writetable(RQ_predictions,['../../models/prot_constrained/' name '/RQ_predictions_' conditions{i} '.txt'],'Delimiter','\t')
%move prot_abundance file back
try
    movefile ../../Databases/prot_abundance_temp.txt ../../Databases/prot_abundance.txt
catch
	disp('prot_abundance.txt file not found in Databases folder') 
end
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
function model = DataConstrains(model,compounds,bounds,flexBounds)
if ~isempty(compounds)
    disp('Constraining byproducts exchange fluxes with fermentation data')
    for i=1:length(compounds)
        %Get exchange rxn index
        if ~strcmpi(compounds{i},'oxygen')
            rxnName = [compounds{i} ' exchange'];
        else
            rxnName = [compounds{i} ' exchange (reversible)'];
        end
        BPindex = find(strcmpi(model.rxnNames,rxnName));
        if ~isempty(BPindex)
            disp([compounds{i} ' exchange has been constrained to: ' num2str(bounds(i)) ' [mmol/gDw h]'])
            %Allow some flexibility 
            model = setParam(model,'ub',BPindex,flexBounds(1)*bounds(i));
            if numel(flexBounds)>1
                model = setParam(model,'lb',BPindex,flexBounds(2)*bounds(i));
            end
        else
            disp(['No exchange rxn for ' compounds{i} ' was found in ecModel'])
        end
    end
end
disp(' ')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function oxPhos = getOxPhosRxnIDs(model,parameters)
if isfield(parameters,'oxPhos')
    oxPhos = [];
    for i=1:length(parameters.oxPhos)
        ID = parameters.oxPhos{i};
        isoEnzymes = model.rxns(find(contains(model.rxns,ID)));
        isoEnzymes = isoEnzymes(~contains(isoEnzymes,'arm_'));
        oxPhos = [oxPhos; isoEnzymes];
    end
else
    oxPhos = [];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeProtCounts(condition,prot_parameters,name)
initial_prots  = length(prot_parameters{1});
filtered_prots = length(prot_parameters{2});
matched_prots  = length(prot_parameters{3});
model_prots    = length(prot_parameters{4});
mass_coverage  = prot_parameters{5};
disp(['The total number of proteins in the dataset is:                ' num2str(initial_prots)])
disp(['The total number of proteins in the filtered dataset is:       ' num2str(filtered_prots)])
disp(['The total number of filtered proteins present in the model is: ' num2str(matched_prots)])
disp(['The mass ratio between measured and unmeasured protein is:     ' num2str(mass_coverage)])

T = table(initial_prots,filtered_prots,matched_prots,model_prots,mass_coverage);
writetable(T,['../../models/prot_constrained/' name '/prot_counts_' condition '.txt'],'Delimiter','\t')
end
