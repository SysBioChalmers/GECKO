%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function FluxRange = MAXmin_Optimizer(model,indexes)
%  
% Get a model and the index(es) of the rxns to maximize and minimize. If
% both optimizations were feasible, then a FV range is returned.
%
% Ivan Domenzain.      Last edited: 2019-02-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FluxRange = MAXmin_Optimizer(model,indexes,bounds,tol)
    FluxRange = [];
    rev       = (length(indexes)>1);
    %Index is a 2 cells array when the rxn is reversible
    for i=1:length(indexes)
        %%% Maximization 
        %Set objective function to maximization of the desired flux
        model.c             = zeros(length(model.c),1);
        model.c(indexes(i)) = 1;
        %Fixes the flux for the backward rxn for irrev models
        Temp_model = setBounds(model,indexes,rev,bounds);
        solution   = solveLP(Temp_model);
        %If Maximization was feasible, then proceed to minimization %%%%
        if ~isempty(solution.f)
            maxFlux             = solution.x(indexes(i));
            model.c             = zeros(length(model.c),1);
            model.c(indexes(i)) = -1;
            %Fixes the flux for the forward rxn for irrev models
            Temp_model = setBounds(model,indexes,false,maxFlux);
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
function model = setBounds(model,indexes,rev,bounds)
    if length(indexes)>1
        if rev 
            %Sets an upper bound for backward reaction to avoid artificially 
            %induced high variability
            model.ub(indexes(2)) = bounds(2);
        else
            %Sets an upper bound for forward reaction with the maximum flux 
            %obtained in the previous step to avoid artificially induced
            %high variability
            model.ub(indexes(1)) = bounds;
        end
    end
end
