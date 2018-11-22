function [FVA_Dists,indexes,stats] = comparativeFVA(model,ecModel,c_source,chemostat,tol,blockedMets)
% comparativeFVA
%  
% This function goes through each of the rxns in a metabolic model and
% gets its flux variability range, then the rxn is mapped into an EC
% version of it to perform the correspondent variability analysis and
% finally compares and plots the cumulative flux variability distributions. 
%
%   model       MATLAB GEM structure
%   ecModel     MATLAB ecGEM structure
%   c_source    rxnName for the main carbon source uptake reaction (which
%               is fixed on the model according to the obtained value from 
%               a batch growth optimization with the ecModel)
%   chemostat   TRUE if chemostat conditions are desired
%   blockedMets metNames of metabolites which secretion should be blocked
%   tol         (Opt) numerical tolerance for a flux and variability range 
%               to be considered as zero
%   
%   FVAdists    cell containing the distributions of variability ranges for
%               the original GEM and ecGEM
%   rangeEC     Distribution of variability ranges for the original ecGEM
%   indexes     Indexes (in the original model) of the reactions for which 
%               a feasible variability range was obtained in both models
%   stats       Some statistics of the variability distributions
% 
% usage: [FVA_Dists,indexes,stats] = comparativeFVA(model,ecModel,c_source,chemostat,tol,blockedMets)
% 
% Ivan Domenzain.      Last edited: 2018-11-21

current  = pwd;
rangeGEM = [];
indexes  = [];
range_EC = [];

if nargin<5
    tol = 0;
end
%Get the index for all the non-objective rxns in the original irrevModel
rxnsIndxs    = find(model.c~=1);
%Set minimal glucose media for ecModel
cd ../../scripts
[ecModel,pos] = changeMedia_batch(ecModel,[c_source ' (reversible)'],'Min');

%Block glucose and oxygen production
cd (current)
if nargin>5
    model   = block_production(model,blockedMets,true);
    ecModel = block_production(ecModel,blockedMets,true);
end

%Gets the optimal value for ecirrevModel and fixes the objective value to
%this for both irrevModels
if chemostat
    gRate   = 0.1;
    gIndex  = find(ecModel.c,1);
    %Fix dilution rate
    [~,~, ecModel] = fixObjective(ecModel,true,gRate);
    %Fix minimal carbon source uptake rate
    ecModel = setParam(ecModel,'obj', pos(1), -1);
    [~,~,ecModel] = fixObjective(ecModel,false);
    %Fix minimal total protein usage
    index   = find(contains(ecModel.rxnNames,'prot_pool'));
    ecModel = setParam(ecModel,'obj', index, -1);
    [~,ecFluxDist,ecModel] = fixObjective(ecModel,false);
else
    %Optimize growth
    [gRate,~, ecModel] = fixObjective(ecModel,true);
    %Fix minimal total protein usage
    index   = find(contains(ecModel.rxnNames,'prot_pool'));
    ecModel = setParam(ecModel,'obj', index, -1);
    [~,ecFluxDist,ecModel] = fixObjective(ecModel,false);
end

%Fix carbon source uptake and growth rates for the ecModel on the original model
carbonUptake = ecFluxDist(pos(1));
disp([c_source ': ' num2str(carbonUptake)])
c_source           = find(strcmpi(model.rxnNames,c_source));
model.lb(c_source) = -carbonUptake;
[~,FluxDist,model] = fixObjective(model,true,gRate);

% Get the variability range for each of non-objective reactions in the
% original irrevModel
for i=1:length(rxnsIndxs) 
    indx        = rxnsIndxs(i);
    rxnID       = model.rxns(indx);
    %mappedIndxs = rxnMapping(rxnID,irrevModel,false);
    FixedValues = [];
    range       = MAXmin_Optimizer(model,indx,FixedValues,tol);
    %If max and min were feasible then the optimization proceeds with
    %the ecModel
    if ~isempty(range)
        if tol>0
            condition = ~(range<tol & abs(FluxDist)<tol);
        else
            condition = true;
        end
        
        if condition 
            mappedIndxs = rxnMapping(rxnID,ecModel,true);
            FixedValues = ecFluxDist(mappedIndxs);
            rangeEC     = MAXmin_Optimizer(ecModel,mappedIndxs,FixedValues,tol);
            if ~isempty(rangeEC)
                rangeGEM = [rangeGEM; range];
                range_EC = [range_EC; rangeEC];
                indexes  = [indexes; indx];
            end
        end
    end
    higher = (rangeEC-range)/range; 
    disp(['ready with #' num2str(i) ', relative variability reduction:' num2str(higher)])
end
%Plot FV cumulative distributions
FVA_Dists  = {rangeGEM, range_EC};
legends    = {'model', 'ecModel'};
titleStr   = 'Flux variability cumulative distribution';
[~, stats] = plotCumDist(FVA_Dists,legends,titleStr,false);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptimalValue, optFluxDist, irrevModel] = fixObjective(irrevModel,fixed,priorValue)
% Optimize and fixes objective value for GEM
objIndx  = find(irrevModel.c~=0);
if fixed 
    factor = 0.99999;
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
end
disp(['The optimal value is ' num2str(OptimalValue)])
end





