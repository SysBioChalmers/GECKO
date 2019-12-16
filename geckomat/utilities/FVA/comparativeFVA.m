function [FVA_Dists,indexes,blocked,stats] = comparativeFVA(model,ecModel,c_source,chemostat,tol,optGRate)
% comparativeFVA
%  
% This function goes through each of the rxns in a metabolic model and
% gets its flux variability range, then the rxn is mapped into an EC
% version of it to perform the correspondent variability analysis and
% finally compares and plots the cumulative flux variability distributions. 
%
%   model       MATLAB GEM structure (reversible model), constrained with 
%               the desired culture medium constraints and biomass
%               pseudorxn as an objective function.
%   ecModel     MATLAB ecGEM structure, constrained with the desired culture
%               medium biomass pseudorxn as an objective function.
%   c_source    rxn ID (in model) for the main carbon source uptake reaction.
%               The rxn ID should not contain the substring "_REV" in order
%               to avoid any confusion when mapping it to the ecModel
%   chemostat   TRUE if chemostat conditions are desired
%   tol         numerical tolerance for a flux and variability range 
%               to be considered as zero
%   optGRate    TRUE if optimal Growth rate should be obtained and fixed
%               for each of the compared models, otherwise the optimal
%               value for the ecModel is also fixed in model. DEFAULT =
%               TRUE
%   
%   FVAdists    cell containing the distributions of variability ranges for
%               the original GEM and ecGEM
%   rangeEC     Distribution of variability ranges for the original ecGEM
%   indexes     Indexes (in the original model) of the reactions for which 
%               a feasible variability range was obtained in both models
%   blocked     rxn indexes for the blocked rxns (cannot carry flux).
%   stats       Some statistics of the variability distributions
% 
% usage: [FVA_Dists,indexes,stats] = comparativeFVA(model,ecModel,c_source,chemostat,tol,blockedMets)
% 
% Ivan Domenzain.      Last edited: 2019-12-16


if nargin<6
    optGRate = true;
    if nargin<5
        tol = 1E-12;
        if nargin<4
            chemostat = false;
        end
    end
end
%Initialize variables
rangeGEM = [];
indexes  = [];
range_EC = [];
%Gets main carbon source uptake reaction index from both models
posCUR_ec = find(strcmpi(ecModel.rxns,[c_source '_REV']));
posCUR    = strcmpi(model.rxns,c_source);
%Get original objective function index for model
posObj    = find(model.c);
%Gets the optimal value for ecirrevModel and fixes the objective value to
%this for both models
if chemostat
    Drate = 0.1;
    %Fix dilution rate
    [gRate,~, ecModel] = fixObjective(ecModel,true,Drate);
    %Fix minimal carbon source uptake rate in both models
    ecModel         = setParam(ecModel,'obj', posCUR_ec, -1);
    [CUR,~,ecModel] = fixObjective(ecModel,true);
    model           = setParam(model,'lb', posCUR, (-1.0001*CUR));
    model           = setParam(model,'ub', posCUR, (-0.9999*CUR));
else
    %Optimize growth
    [gRate,~, ecModel] = fixObjective(ecModel,true);
    if optGRate
        [gRate,~,~] = fixObjective(model,true);
    end
end
%Get a parsimonious flux distribution for the ecModel (minimization of
%total protein usage)
Pool_index         = find(contains(ecModel.rxnNames,'prot_pool'));
ecModel            = setParam(ecModel,'obj', Pool_index, -1);
[~,ecFluxDist,~]   = fixObjective(ecModel,false);
model              = setParam(model,'obj', posObj, 1);
[~,FluxDist,model] = fixObjective(model,true,gRate);
%Get the index for all the reactions that can carry a flux in the original
%model and then run FVA on that subset
disp('Identifying reactions that can carry a non-zero flux')
rxnsIndxs = haveFlux(model,tol);
rxnsIndxs = find(model.c~=1 & rxnsIndxs);
blocked   = rxnsIndxs(rxnsIndxs==0);
%Get the variability range for each of the flux carrier reactions
if ~isempty(FluxDist) && ~isempty(rxnsIndxs)
    disp('Performing Flux Variability Analysis')
    for i=1:length(rxnsIndxs)
        indx  = rxnsIndxs(i);
        rxnID = model.rxns(indx);
        range = MAXmin_Optimizer(model,indx,1000,tol);
        rev   = false;
        if model.rev(indx) ==1
            rev = true;
        end
        %If max and min were feasible then the optimization proceeds with
        %the ecModel
        if ~isempty(range)
            %Get the correspondent index(es) for the i-th reaction in the
            %ecModel
            mappedIndxs = rxnMapping(rxnID,ecModel,rev);
            %Get bounds from the optimal distribution to avoid artificially
            %induced variability
            bounds  = ecFluxDist(mappedIndxs);
            rangeEC = MAXmin_Optimizer(ecModel,mappedIndxs,bounds,0);
            if ~isempty(rangeEC)
                rangeGEM = [rangeGEM; range];
                range_EC = [range_EC; rangeEC];
                indexes  = [indexes; indx];
                disp(['ready with #' num2str(i) ' | model Variability: ' num2str(range) ' | ecModel variability: ' num2str(rangeEC)])
            end
        end
    end
else
    warning('The metabolic model is unfeasible under the provided constraints')
end
%Plot FV cumulative distributions
FVA_Dists  = {rangeGEM, range_EC};
legends    = {'model', 'ecModel'};
titleStr   = 'Flux variability cumulative distribution';
[~, stats] = plotCumDist(FVA_Dists,legends,titleStr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptimalValue, optFluxDist, model] = fixObjective(model,fixed,priorValue)
% Optimize and fixes objective value for GEM
objIndx  = find(model.c~=0);
if nargin < 3
    sol = solveLP(model);
    priorValue = sol.x(objIndx);
end

if fixed 
    factor = 0.9999;
else
    factor = 0;
end

model.lb(objIndx) = factor*priorValue;
model.ub(objIndx) = 1.0001*priorValue;

sol = solveLP(model);
if ~isempty(sol.f)
    OptimalValue = sol.x(objIndx);
    optFluxDist  = sol.x;
end
disp(['The optimal value for ' model.rxns{objIndx} ' is ' num2str(OptimalValue)])
end