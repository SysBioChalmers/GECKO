function generate_protModels(ecModel,grouping,flexFactor,oxPhosIDs,ecModel_batch)
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
%   oxPhosIDs     (cell) IDs for the oxidative phosphorylation related. The
%                 mean expression level across subunits is taken as the 
%                 abundance levelapplied to each of the enzyme complexes 
%                 in these reactions to avoid overconstraining.
%   ecModel_batch (structure) ecModel MATLAB structure with total protein pool
%                 constraint, if provided the process of fitting GAM for
%                 each new protein content is speed-up.
%
% Usage: generate_protModels(ecModel,conditions,grouping,byProducts)
%
% Last modified.  Ivan Domenzain 2019-09-10

close all
current = pwd;

%This funcion allows for flexibilization of protein absolute abundances in 
%case that ecModelP is not feasible using the automatically flexibilized 
%data, if flex factor is not specified then a factor of 1 is assumed.
if nargin<5
    ecModel_batch = [];
    if nargin<4
        oxPhosIDs = [];
        if nargin<3
            flexFactor = 1.05;
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

%For each condition create a protein constrained model
for i=1:length(conditions)
    cd (current)
    disp(conditions{i})
    %Extract data for the i-th condition
    abundances = cell2mat(absValues(1:grouping(i)));
    absValues  = absValues(grouping(i)+1:end);
    %Filter data
    [pIDs, abundances] = filter_ProtData(uniprotIDs,abundances,1.96,true);
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
            [~,GAM]   = scaleBioMass(tempModel,Ptot(i),[],true);
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
    f              = 0.5; %Protein mass in model/Total theoretical proteome
    flexGUR        = flexFactor*GUR(i);
    disp(['Incorporation of proteomics constraints for ' conditions{i} ' condition'])
    [ecModelP,~,~] = constrainEnzymes(ecModelP,f,GAM,Ptot(i),pIDs,abundances,Drate(i),flexGUR);
    
    %Set chemostat conditions constraints and fit NGAM
    cd (current)
    %NGAM interval for fitting
    interval = [0 5];
    ecModelP = setStressConditions(ecModelP,Drate(i),positionsEC,expData,NGAM,interval);
    
    %Get optimal flux distribution and display exchange fluxes
    solution = solveLP(ecModelP,1);
    if ~isempty(solution.f)
        printFluxes(ecModelP,solution.x)
    end
    %Fix experimental Glucose uptake rate and save models
    cd ..
    ecModelP = setChemostatConstraints(ecModelP,positionsEC,Drate(i),true,0.001,GUR(i));
    save(['../../models/ecModel_Prot_' conditions{i} '.mat'],'ecModelP')
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
condModel          = setChemostatConstraints(model,pos,Drate,minProt);
condModel.c(:)     = 0;
prots              = find(contains(condModel.rxnNames,'prot_'));
protPool           = find(contains(condModel.rxnNames,'prot_pool'));
prots              = prots(prots~=protPool);
condModel.c(prots) = -1;
cd integrate_proteomics
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
            disp([byProducts{i} ' exchange has been constrained with an UB: ' num2str(bounds(i)) ' [mmol/gDw h]'])
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

