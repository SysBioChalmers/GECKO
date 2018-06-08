% function model = flexibilizeProteins(model,gRate)
%
% Takes an ecModel_batch structure and performs a growth rate optimization
% on minimal glucose media, if the model is overconstrained the algorithm
% will iterate flexibilizing the upper bounds for the protein exchange
% reaction until it grows at the specified experimental value
%
% Created.  Ivan Domenzain 2018-05-02
%
function [model,enzUsages,modifications] = flexibilizeProteins(model,gRate,glucUptakeExp)
current    = pwd;
flexFactor = 1000;
flexProts  = {};
% set minimal glucose medium
Csource      = 'D-glucose exchange (reversible)';
glucUptkIndx = find(strcmpi(model.rxnNames,Csource));
cd ../Kcat_sensitivity_analysis
[model,~] = changeMedia_batch(model,Csource,'Min');
%constrain glucose uptake if an experimental measurement is provided
if nargin>2
    model.ub(glucUptkIndx) = 1.01*glucUptakeExp;
end
cd (current)
% get measured protein exchange rxns indexes
measuredIndxs = getMeasuredProtsIndexes(model);
abundances    = model.ub(measuredIndxs);
objIndex      = find(model.c==1);
sol           = solveLP(model,1);
growth        = sol.x(objIndex);
% iterate while growth is underpredicted
while growth<0.99*gRate
    limIndex = findLimitingUBs(model,measuredIndxs,objIndex,flexFactor);
    %Flexibilize the top growth limiting protein on the original eModel
    flexProts          = [flexProts; limIndex];
    model.ub(limIndex) = Inf;%model.ub(limIndex)*flexFactor;
    sol                = solveLP(model);
    if ~isempty(sol.x)
       growth = sol.x(objIndex);
       disp(['Modified ub for: ' model.rxnNames{limIndex} ' gRate: ' num2str(growth)])
    end
end
[model,enzUsages]  = getNewBounds(model,gRate,measuredIndxs,flexProts,objIndex);
modifiedAbundances = model.ub(measuredIndxs);
exchangedProteins  = model.rxnNames(measuredIndxs);
modifications      = getDifferences(abundances,modifiedAbundances,exchangedProteins);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffTable = getDifferences(originalBounds,newBounds,exchangedProteins)
protein_IDs     = {};
previous_values = {};
modified_values = {};
for i=1:length(originalBounds)
    if newBounds(i)>originalBounds(i)
        proteinID       = exchangedProteins{i};
        proteinID       = proteinID(1:((strfind(proteinID,'_exchange'))-1));
        protein_IDs     = [protein_IDs; proteinID];
        previous_values = [previous_values; originalBounds(i)];
        modified_values = [modified_values; newBounds(i)];
    end
end
diffTable = table(protein_IDs,previous_values,modified_values);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxIndex = findLimitingUBs(model,protIndxs,objIndex,flexFactor)
coefficients = zeros(length(protIndxs),1);
solution     = solveLP(model,1);
    % flexibilize ub for every protein exchange in a temporal model
    for i=1:length(protIndxs)
        index               = protIndxs(i);
        tempModel           = model;
        tempModel.ub(index) = tempModel.ub(index)*flexFactor;
        %get a new solution with the flexibilized ub and calculate the
        %effect on the growth rate prediction
        newSol = solveLP(tempModel,1);
        if ~isempty(newSol.f)
            ControlCoeff    = (newSol.x(objIndex)-solution.x(objIndex))/solution.x(objIndex);
            coefficients(i) = abs(ControlCoeff);
        end
    end
    [~,maxIndex] = max(coefficients);
    maxIndex     = protIndxs(maxIndex);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measuredIndxs = getMeasuredProtsIndexes(model)
measuredIndxs  = find(contains(model.rxnNames,'prot_'));
exchange_prots = find(contains(model.rxnNames(measuredIndxs),'_exchange'));
measuredIndxs  = measuredIndxs(exchange_prots(1:end-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,enzUsages] = getNewBounds(model,gRate,protIndxs,flexProts,gPos)
%Now that the model is growing at least at the specified dilution rate
%lets fix the growth rate and minimize enzymes usage
model.lb(gPos) = 0.999*gRate;
model.ub(gPos) = 1.001*gRate;
model.c(:)           = 0;
protIndexes          = find(contains(model.rxnNames,'prot_'));
model.c(protIndexes) = -1;
optSolution          = solveLP(model,1);
optSolution          = optSolution.x;
enzUsages            = zeros(length(protIndxs),1);
for i=1:length(protIndxs)
    index = protIndxs(i);
    %If protein was flexibilized set its upper bound to the simulated
    %concentration
    if ismember(protIndxs(i),flexProts) && optSolution(index)>0
        model.ub(index) = optSolution(index);
        %newSol          = solveLP(model,1);
    end
    enzUsages(i) = optSolution(index)/model.ub(index);
end
end
    
    
