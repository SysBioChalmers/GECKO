function [FVA_Dists,indexes,blocked,stats] = comparativeFVA(model,ecModel,c_source,chemostat,tol)
% comparativeFVA
%  
% This function goes through each of the rxns in a metabolic model and
% gets its flux variability range, then the rxn is mapped into the EC
% version of the model to perform the correspondent variability analysis and
% finally compares and plots the cumulative flux variability distributions. 
%
%   model       MATLAB GEM structure (reversible model), constrained with 
%               the desired culture medium constrains
%   ecModel     MATLAB ecGEM structure, constrained with the desired culture
%               medium
%   c_source    rxn ID (in model) for the main carbon source uptake reaction.
%               The rxn ID should not contain the substring "_REV" in order
%               to avoid any confusion when mapping it to the ecModel
%   chemostat   TRUE if chemostat conditions are desired
%   tol         numerical tolerance for a flux and variability range 
%               to be considered as zero (default 1E-12) 
%   FVAdists    cell containing the distributions of variability ranges for
%               the original GEM and ecGEM
%   rangeEC     Distribution of variability ranges for the original ecGEM
%   indexes     Indexes (in the original model) of the reactions for which 
%               a feasible variability range was obtained in both models
%   stats       Some statistics of the variability distributions
% 
% usage: [FVA_Dists,indexes,stats] = comparativeFVA(model,ecModel,c_source,chemostat,tol,blockedMets)
% 
% Ivan Domenzain.      Last edited: 2019-05-13

if nargin<5
    tol = 1E-12;
end
rangeGEM = [];
indexes  = [];
blocked  = [];
range_EC = [];
%Get the index for all the non-objective rxns in the original irrevModel
rxnsIndxs = find(model.c~=1);
%Constraint all rxns in ecModel to the positive domain
ecModel.lb = zeros(length(ecModel.lb),1);
%Gets main carbon source uptake reaction index from both models
posCS    = strcmpi(model.rxns,c_source);
posCS_ec = find(strcmpi(ecModel.rxns,[c_source '_REV']));

%If chemostat conditions are desired then a dilution rate is fixed for the
%ecModel, CS uptake rate set as objective for minimization, then its
%optimal value is fixed
if chemostat
    gRate   = 0.1;
    %Fix dilution rate in ecModel
    [~,~, ecModel] = constrainIrrevModel(ecModel,true,gRate);
    %Fix minimal carbon source uptake rate in ecModel
    ecModel        = setParam(ecModel,'obj', posCS_ec, -1);
    [~,~,ecModel]  = constrainIrrevModel(ecModel,false);
else
    %For batch conditions growth is maximized and then fixed 
    [gRate,~, ecModel] = constrainIrrevModel(ecModel,true);
end
%Set minimal total protein usage in ecModel as objective function
index   = find(contains(ecModel.rxnNames,'prot_pool'));
ecModel = setParam(ecModel,'obj', index, -1);
%Fix optimal protein usage and get an optimal flux distribution for ecModel
[~,ecFluxDist,ecModel] = constrainIrrevModel(ecModel,true);
%Get optimal Csource uptake rate from optimal ecModel flux distribution
Cuptake  = ecFluxDist(posCS_ec);
disp([c_source ': ' num2str(Cuptake)])
%Constrain model with values from the optimal ecModel flux distribution
[model,FluxDist] = constrainModel(model,gRate,posCS,chemostat,Cuptake);
%Get the variability range for each of the non-objective reactions in the
%original model
if ~isempty(FluxDist) & ~isempty(rxnsIndxs)
    for i=1:length(rxnsIndxs)
        indx    = rxnsIndxs(i);
        fluxVal = FluxDist(indx);
        rxnID   = model.rxns(indx);
        rev     = false;
        if model.rev(indx) ==1
            rev = true;
        end
        bounds = [];
        range  = MAXmin_Optimizer(model,indx,bounds,tol);
        %If max and min were feasible then the optimization proceeds with
        %the ecModel
        if ~isempty(range)
            %MAX-min proceeds for the ecModel if the FV range and optimal flux
            %value are non-zero for the original model
            if ~(range<tol & abs(fluxVal)<tol)
                %Get the correspondent index(es) for the i-th reaction in the
                %ecModel
                mappedIndxs = rxnMapping(rxnID,ecModel,rev);
                %Get bounds from the optimal distribution to avoid artificially
                %induced variability
                bounds      = ecFluxDist(mappedIndxs);
                rangeEC     = MAXmin_Optimizer(ecModel,mappedIndxs,bounds,tol);
                if ~isempty(rangeEC)
                    rangeGEM = [rangeGEM; range];
                    range_EC = [range_EC; rangeEC];
                    indexes  = [indexes; indx];
                    %disp(['ready with #' num2str(i) ' // model Variability: ' num2str(range) ' // ecModel variability: ' num2str(rangeEC)])
                end
            else
                blocked  = [blocked; indx];
            end
        end
    end
else
    warning('The original model is unfeasible under the provided constraints')
end
%Plot FV cumulative distributions
FVA_Dists  = {rangeGEM, range_EC};
legends    = {'model', 'ecModel'};
titleStr   = 'Flux variability cumulative distribution';
[~, stats] = plotCumDist(FVA_Dists,legends,titleStr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,FluxDist] = constrainModel(model,gRate,posCS,chemostat,Cuptake)
%Fix optimal carbon source uptake and growth rate values from the ecModel in the 
% original model (Convention: GEMs represent uptake fluxes as negative values)
[~,FluxDist,model] = constrainIrrevModel(model,true,gRate);
if chemostat 
    %Set glucose uptake rate as objective to minimize
    model = setParam(model,'obj', posCS, 1);
else
    %Fix glucose uptake rate
    model = setParam(model,'lb', posCS,-1.0001*Cuptake);
    model = setParam(model,'ub', posCS,-0.9999*Cuptake);
end
FluxDist = solveLP(model,1);
FluxDist = FluxDist.x;
if ~isempty(FluxDist)
    Cuptake = FluxDist(posCS);
    model   = setParam(model,'lb', posCS, 1.0001*Cuptake);
    disp(['The optimal value for ' model.rxnNames{posCS} ' is: ' num2str(Cuptake)])
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptimalValue, optFluxDist, irrevModel] = constrainIrrevModel(irrevModel,fixed,priorValue)
% Optimize and fixes objective value for ecModel
objIndx      = find(irrevModel.c~=0);
optFluxDist  = [];
OptimalValue = 0;
if fixed 
    factor = 0.9999;
else
    factor = 0;
end

if nargin == 3
    irrevModel.lb(objIndx) = factor*priorValue;
    irrevModel.ub(objIndx) = priorValue;
else
    sol = solveLP(irrevModel);
    irrevModel.lb(objIndx) = factor*sol.x(objIndx);
    irrevModel.ub(objIndx) = sol.x(objIndx);
end
sol = solveLP(irrevModel);
if ~isempty(sol.f)
    OptimalValue = sol.x(objIndx);
    optFluxDist  = sol.x;
    disp(['The optimal value for ' irrevModel.rxnNames{objIndx} ' is: ' num2str(OptimalValue)])
else
    disp('Non feasible simulation')
end
end