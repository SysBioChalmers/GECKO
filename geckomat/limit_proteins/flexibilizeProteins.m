function [model,enzUsages,modifications,gRate] = flexibilizeProteins(model,gRate,c_UptakeExp,c_source)
% flexibilizeProteins
%   Function that takes an ecModel with proteomic constraints and, if it is
%   overconstrained with respect to the provided experimental growth rate,
%   iterates finding the top growth-limiting enzyme or enzyme complex. The
%   exchange rate upper bound for each of the identified enzyme or subunits 
%   is then set to infinity and after all iterations they are set equal to 
%   the usage value provided by an overall enzyme usage minimization 
%   simulation (subject to the provided growth rate and nutrient uptake 
%   constraints).
%
%   model           ecModel with proteomic constraints (individual enzyme
%                   levels)
%   gRate           Minimum growth rate the model should grow at [1/h]. For
%                   finding the growth reaction, GECKO will choose the
%                   non-zero coeff in the objective function.
%   c_UptakeExp     (Opt) Experimentally measured glucose uptake rate 
%                   [mmol/gDw h]
%	c_source        (Opt) The name of the exchange reaction that supplies
%                   the model with carbon.
%
%   model           ecModel with calibrated enzyme usage upper bounds
%   enzUsages       Calculated enzyme usages after final calibration 
%                   (enzyme_i demand/enzyme_i upper bound)
%   modifications   Table with all the modified values 
%                   (Protein ID/old value/Flexibilized value)
%
%   Usage: [model,enzUsages,modifications] = flexibilizeProteins(model,gRate,c_UptakeExp,c_source)
%
%   Benjamin J. Sanchez     2018-12-11
%   Ivan Domenzain          2020-02-06
%

flexFactor    = 100;
flexProts     = {};
enzUsages     = [];
modifications = {};
PoolIndex     = find(contains(model.rxnNames,'prot_pool_exchange'));
%constrain glucose uptake if an experimental measurement is provided
if nargin > 2
    glucUptkIndx           = strcmp(model.rxnNames,c_source);
    model.ub(glucUptkIndx) = c_UptakeExp;
end
% get measured protein exchange rxns indexes
measuredIndxs = getMeasuredProtsIndexes(model);
if ~isempty(measuredIndxs)
    abundances = model.ub(measuredIndxs);
    objIndex   = find(model.c==1);
    sol        = solveLP(model,1);
    growth     = sol.x(objIndex);
    difference = gRate-growth;
    % iterate while growth is underpredicted
    while growth<gRate && difference>(-1E-6)
        tempModel = model;
        [limIndex,flag] = findLimitingUBs(model,measuredIndxs,flexFactor,1);
        if ~flag
            [limIndex,~] = findLimitingUBs(model,measuredIndxs,flexFactor,2);
        end
        if ~isempty(limIndex)
            %Flexibilize the top growth limiting protein on the original ecModel
            flexProts = [flexProts; model.rxns(limIndex)];
            tempModel.ub(limIndex) = Inf;
            sol = solveLP(tempModel);
            if ~isempty(sol.x)
            	growth = sol.x(objIndex);
                for j=1:length(limIndex)
                    indx = limIndex(j);
                    disp(['Modified ub for: ' model.rxnNames{indx} ' gRate: ' num2str(growth)])
                end
                %Update model structure
                model = tempModel;
            else
               %In case that the resulting model is a non-functional one
               %then proceed with a suboptimal growth rate (this makes the 
               %while loop to break)
               warning(['Unfeasible flexibilization of ' model.rxnNames{indx} ' UB'])
               gRate = growth; 
            end
        else
            %In case that no limiting enzymes have been found then proceed 
            %with a suboptimal growth rate (this makes the while loop to break)
            warning('No limiting enzymes were found')
            gRate = growth;
        end
    end
    [model,enzUsages]  = getNewBounds(model,gRate,measuredIndxs,flexProts,objIndex);
    modifiedAbundances = model.ub(measuredIndxs);
    exchangedProteins  = model.rxnNames(measuredIndxs);
    modifications      = getDifferences(abundances,modifiedAbundances,exchangedProteins,model);
    %Get total flexibilized protein mass
    totalFlexMass = sum(modifications.flex_Mass);
    %Allow 100% of saturation for protein pool enzymes
    poolIndex = contains(model.rxnNames,'prot_pool_exchange');
    Ppool     = model.ub(poolIndex);   
    tempModel = setParam(model,'ub',poolIndex,2*Ppool);
    %Get the minimal pool usage
    tempModel = setParam(tempModel,'obj',poolIndex,-1);
    tempModel = setParam(tempModel,'lb',objIndex,0.99*gRate);
    sol = solveLP(tempModel);
    if ~isempty(sol.x)
        %Compare minimal pool usage with flex prot mass
        minPool = sol.x(poolIndex);
        if (Ppool - minPool)>=0.5*totalFlexMass
            tempModel.ub(poolIndex) = Ppool-0.5*totalFlexMass;
            disp('Successful protein flexibilization')
            model = tempModel;
        else
            warning('Flexibilized protein mass exceeds available protein pool')
        end
    else
        warning('Unfeasible protein flexibilization')
    end        
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffTable = getDifferences(originalBounds,newBounds,exchangedProteins,model)
protein_IDs     = {};
previous_values = [];
modified_values = [];
Mweights        = [];
for i=1:length(originalBounds)
    if newBounds(i)>originalBounds(i)
        proteinID   = exchangedProteins{i};
        proteinID   = proteinID(1:((strfind(proteinID,'_exchange'))-1));
        protein_IDs = [protein_IDs; proteinID];

        %Get protein MW
        index = find(strcmpi(model.enzymes,strrep(proteinID,'prot_','')));
        if ~isempty(index)
            Mweights = [Mweights; model.MWs(index)];
        end
        %Get previous and modified usage values [mmol/gDw h]
        previous_values = [previous_values; originalBounds(i)];
        if newBounds(i)~=Inf
            modified_values = [modified_values; newBounds(i)];
        else
            modified_values = [modified_values; originalBounds(i)];
        end
    end
end
%Calculate flexibilized mass for each protein [mmol/gDw h]*[g/mmol]>g/gDwh
flex_Mass = (modified_values-previous_values).*Mweights;
diffTable = table(protein_IDs,previous_values,modified_values,flex_Mass);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measuredIndxs = getMeasuredProtsIndexes(model)
measuredIndxs  = intersect(find(contains(model.rxnNames,'prot_')),find(contains(model.rxnNames,'_exchange')));
measuredIndxs  = measuredIndxs(1:end-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,usagesTable] = getNewBounds(model,gRate,protIndxs,flexProts,gPos)
%Now that the model is growing at least at the specified dilution rate
%lets fix the growth rate and minimize enzymes usage
objectiveVector = model.c;
model.lb(gPos)  = 0.9999*gRate;
model.c(:)      = 0;
protNames       = model.rxnNames(protIndxs);
pool_Index      = contains(model.rxnNames,'prot_pool_');
%Forces flux to split over measured enzymes
model.c(pool_Index)  = -1;
optSolution          = solveLP(model,1);
optSolution          = optSolution.x;
enzUsages            = zeros(length(protIndxs),1);
for i=1:length(protIndxs)
    index = protIndxs(i);
    name  = model.rxns(index);
    %If protein was flexibilized set its upper bound to the simulated
    %concentration
    if ismember(name,flexProts) && optSolution(index)>0
        model.ub(index) = optSolution(index);
    end
    enzUsages(i) = optSolution(index)/model.ub(index);
end
model.c     = objectiveVector;
usagesTable = table(protNames,enzUsages,'VariableNames',{'prot_IDs' 'usage'});
end
    
    
