%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function FluxRange = MAXmin_Optimizer(model,indexes)
%  
% Get a model and the index(es) of the rxns to maximize and minimize. If
% both optimizations were feasible, then a FV range is returned.
%
% Ivan Domenzain.      Last edited: 2018-03-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FluxRange = MAXmin_Optimizer(model,indexes,FixedValues,tol)
    FluxRange = [];
    %Index is a 2 cells array when the rxn is a splitted reversible rxn
    for i=1:length(indexes)
        %%% Maximization 
        %Set objective function to maximization of the desired flux
        model.c             = zeros(length(model.c),1);
        model.c(indexes(i)) = 1;
        %Fixes the flux for the backward rxn for irrev models
        Temp_model = fixFluxValue(model,indexes,1,FixedValues);
        solution   = solveLP(Temp_model);
        %If Maximization was feasible, then proceed to minimization %%%%
        if ~isempty(solution.f)
            maxFlux             = solution.x(indexes(i));
            model.c             = zeros(length(model.c),1);
            model.c(indexes(1)) = -1;
            %Fixes the flux for the forward rxn for irrev models
            Temp_model = fixFluxValue(model,indexes,-1,maxFlux);
            solution   = solveLP(Temp_model);
            if ~isempty(solution.f)
                minFlux   = solution.x(indexes(i));
                FluxRange = [FluxRange; maxFlux-minFlux];
            end   
        end
    end
    
    if ~isempty(FluxRange)
        FluxRange = sum(FluxRange);
        if FluxRange<tol
            FluxRange = 0;
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = fixFluxValue(model,indexes,coeff,FixedValues)
    if length(indexes)>1
        if coeff == 1
            %Sets an upper bound for backward reaction to avoid artificially 
            %high variability
            model.ub(indexes(2)) = FixedValues(2);
        else
            %Sets an upper bound for forward reaction with the maximum flux 
            %obtained in the previous step to avoid artificially high variability
            model.ub(indexes(1)) = FixedValues;
        end
    end
end
