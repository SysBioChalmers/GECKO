function generate_protModels(ecModel,grouping,name,ecModel_batch,parameters)
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
%   ecModel_batch (structure) ecModel MATLAB structure with total protein pool
%                 constraint.
%
% Usage:  generate_protModels(ecModel,grouping,name,ecModel_batch)
%
% Last modified.  Ivan Domenzain 2020-04-06
close all
current = pwd;

% Save geckomat path as a parameter to return after use organism specific
% scripts
cd ../..
parameters.geckomat = pwd;

if nargin < 4
   error('ecModel_batch must be specified');
end

if nargin < 5
    error('Not enough input arguments, parameters are required');
end

%This funcion allows for flexibilization of protein absolute abundances in
%case that ecModelP is not feasible using the automatically flexibilized 
%data, if flex factor is not specified then a factor of 1 is assumed.

Ptot_model = parameters.Ptot;
c_source   = parameters.c_source;
growthRxn  = parameters.exch_names{1};
NGAM       = parameters.NGAM;
GAM        = [];

% Validtae if there is a defined path where organism/model specific scripts are
if ~isfield(parameters, 'customPath')
    [~,values] = fileattrib('../../custom/');
    % Add to parameters the full path
    parameters.customPath = [values.Name '/' name '/'];
else
    [~,values] = fileattrib(parameters.customPath);
    %Update to a full path
    parameters.customPath = [values.Name '/' name '/'];
end

% Validtae if there is a defined path where save the output
if ~isfield(parameters, 'outputPath')
    [~,values] = fileattrib('../../ecModels/prot_constrained/');
    % Add to parameters the full path
    parameters.outputPath = [values.Name '/'  name '/'];
else
    [~,values] = fileattrib(parameters.outputPath);
    %Update to a full path
    parameters.outputPath = [values.Name '/prot_constrained/' name '/'];
end

%Flexibilization factor for carbon source uptake rate (needed for
%flexibilizeProteins step in constrainEnzymes).
flexFactor = 1.05;

%Get oxPhos related rxn IDs
oxPhos = getOxPhosRxnIDs(ecModel,parameters);

%Get rxn IDs for carbon source uptake and growth reaction. Since positions
%can change. it is better to used the rxnID.
c_source_ID = ecModel.rxns(find(strcmpi(ecModel.rxnNames,c_source)));
growthRxn_ID = ecModel.rxns(find(strcmpi(ecModel.rxnNames,growthRxn)));

%Remove prot_abundance.txt  and relative_proteomics.txt files
%(for f factor calculation)
try
    movefile([parameters.customPath 'data/prot_abundance.txt'], ...
        [parameters.customPath 'data/prot_abundance_temp.txt']);
catch
    disp('prot_abundance.txt file not found in Databases folder') 
end
%Load absolute proteomics dataset [mmol/gDw]
%and fermentation data (GUR, OUR, CO2 production, byProducts, Ptot, Drate)
cd(current)
[uniprotIDs,absValues,fermData,byProducts] = load_Prot_Ferm_Data(grouping,parameters);
conditions = fermData.conds;
Ptot       = fermData.Ptot;
Drate      = fermData.Drate;
GUR        = fermData.GUR;
CO2prod    = fermData.CO2prod;
OxyUptake  = fermData.OxyUptake;
byP_flux   = fermData.byP_flux;
%For each condition create a protein constrained model
for i=1:length(conditions)
    cd(current)
    %create specific subfolder for ecModelProts output files if not present already
    if ~isfolder([parameters.outputPath conditions{i}])
        mkdir([parameters.outputPath conditions{i}])
    end

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
    cd([parameters.customPath 'scripts'])
    ecModelP  = changeMedia_batch(ecModel,c_source);
    tempModel = changeMedia_batch(ecModel_batch,c_source);

    %If the relative difference between the ecModel's protein content and
    %the Ptot for i-th condition is higher than 5% then biomass should be
    %rescaled and GAM refitted to this condition.
    %For fitting GAM a functional model is needed therefore an ecModel with
    %total protein pool constraint should be used
    Prot_diff = abs(Ptot_model-Ptot(i))/Ptot_model;
    if Prot_diff>=0.05
        [~,GAM] = scaleBioMass(tempModel,Ptot(i),parameters,[],true);
        %Then the GAM and new biomass composition are set in ecModelP, which
        %is not functional yet but should be used for incorporation of
        %proteomics data
        cd([parameters.customPath 'scripts'])
        ecModelP = scaleBioMass(ecModelP,Ptot(i),parameters,GAM,true);
        disp(' ')
    end
    %Block production of non-observed metabolites before data incorporation
    %and flexibilization
    expData  = [GUR(i),CO2prod(i),OxyUptake(i)];
    flexGUR  = flexFactor*GUR(i);
    ecModelP = DataConstrains(ecModelP,byProducts,byP_flux(i,:),1.1);
    %Get a temporary model structure with the same constraints to be used
    %for minimal enzyme requirements analysis. For all measured enzymes
    %(present in the dataset) a minimal usage is obtained from a FBA
    %simulation with the ecModel_batch, then 
    tempModel       = DataConstrains(tempModel,byProducts,byP_flux(i,:),1.1);
    tempModel       = setParam(tempModel,'ub',c_source_ID,flexGUR);
    [matchedEnz,iA] = intersect(pIDs,tempModel.enzymes);
    enzModel        = setParam(tempModel,'lb',growthRxn_ID,Drate(i));
    for j=1:length(matchedEnz)
        rxnIndex  = find(contains(tempModel.rxnNames,matchedEnz{j}));
        tempModel = setParam(enzModel,'obj',rxnIndex,-1);
        tempSol   = solveLP(tempModel);
        %Compare enzyme minimum usage with abundance value
        if (tempSol.x(rxnIndex)-abundances(iA(j)))>0
            %Flexibilize limiting values
            disp(['Limiting abundance found for: ' matchedEnz{j} '/Previous value: ' num2str(abundances(iA(j))) ' /New value: ' num2str(tempSol.x(rxnIndex))])
            abundances(iA(j)) = 1.01*tempSol.x(rxnIndex);
        end
        enzIndex = find(contains(tempModel.enzymes,matchedEnz{j}));
    end
    %Get model with proteomics
    f       = 1; %Protein mass in model/Total theoretical proteome
    disp(['Incorporation of proteomics constraints for ' conditions{i} ' condition'])
    cd([parameters.geckomat '/limit_proteins'])
    [ecModelP,usagesT,modificationsT,~,coverage] = constrainEnzymes(ecModelP,parameters,f,GAM,Ptot(i),pIDs,abundances,Drate(i),flexGUR);
    matchedProteins = usagesT.prot_IDs;
    prot_input = {initialProts filteredProts matchedProteins ecModel.enzymes coverage};
    writeProtCounts(prot_input,parameters.outputPath,conditions{i});
    cd(current)
    %NGAM interval for fitting
    interval = [0 5];
    %Get indexes for carbon source uptake and growth reaction
    positionsEC(1) = find(strcmpi(ecModelP.rxnNames,c_source));
    positionsEC(2) = find(strcmpi(ecModelP.rxnNames,growthRxn));
    ecModelP = setStressConditions(ecModelP,parameters,Drate(i),positionsEC,expData,NGAM,interval);
    %Fix experimental Glucose uptake rate and save models
    cd ..
    ecModelP = setChemostatConstraints(ecModelP,positionsEC,Drate(i),true,0.01,GUR(i));
    %Get optimal flux distribution and display exchange fluxes
    solution = solveLP(ecModelP,1);
    if ~isempty(solution.f)
        fileFluxes = [parameters.outputPath conditions{i} '/fluxes_Exch.txt'];
        printFluxes(ecModelP,solution.x,true,1E-4,fileFluxes)
    end
    cd ../change_model
    [~,~,version] = preprocessModel(ecModelP,'','');
    cd(parameters.geckomat)
    saveECmodel(ecModelP,'RAVEN',name,version,parameters.outputPath,conditions{i});
    %save .txt file
    writetable(usagesT,[parameters.outputPath conditions{i} '/enzymeUsages.txt'],'Delimiter','\t')
    writetable(modificationsT,[parameters.outputPath conditions{i} '/modifiedEnzymes.txt'],'Delimiter','\t')
end

%Remove prot_abundance.txt  and relative_proteomics.txt files
%(for f factor calculation)
try
    movefile([parameters.customPath 'data/prot_abundance_temp.txt'], ...
        [parameters.customPath 'data/prot_abundance.txt']);
catch
    disp('prot_abundance_temp.txt file not found in Databases folder') 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function condModel = setStressConditions(model,parameters,Drate,pos,expData,NGAM,interval)
minProt = true;
ecFlag  = true;
%Rescale biomass set fitted GAM for std conditions and fit non-growth
%associated maintenance for stress conditions
cd ..
condModel = setChemostatConstraints(model,pos,Drate,minProt,0.01);
cd integrate_proteomics
disp(' ')
condModel = fitNGAM(condModel,parameters,NGAM,expData,interval,ecFlag);
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
function writeProtCounts(prot_parameters,path,condition)
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
writetable(T,[path condition '/prot_counts.txt'],'Delimiter','\t')
end
