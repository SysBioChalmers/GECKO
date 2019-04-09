%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function FluxRange = MAXmin_Optimizer(model,indexes)
%  
% Get a model and the index(es) of the rxns to maximize and minimize. If
% both optimizations were feasible, then a FV range is returned.
%
% Ivan Domenzain.      Last edited: 2019-04-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FluxRange = MAXmin_Optimizer(model,indexes,bounds,tol)
    FluxRange = [];
    %Index is a 2 cells array when the rxn is reversible
    for i=1:length(indexes)
        forward  = (i==1);
        %%% Maximization 
        Temp_model = model;
        %Set objective function to maximization of the desired flux
        Temp_model.c             = zeros(length(Temp_model.c),1);
        Temp_model.c(indexes(i)) = 1;
        %Fixes the flux for the backward rxn for irrev models
        Temp_model = setBounds(Temp_model,indexes,forward,bounds);
        solution   = solveLP(Temp_model);
        %If Maximization was feasible, then proceed to minimization %%%%
        if ~isempty(solution.f)
            maxFlux                  = solution.x(indexes(i));
            Temp_model               = model;
            Temp_model.c             = zeros(length(Temp_model.c),1);
            Temp_model.c(indexes(i)) = -1;
            %Get solution vector for minimization of i-th reaction
            solution                 = solveLP(Temp_model);
            if ~isempty(solution.f)
                minFlux = solution.x(indexes(i));
                range   = abs(maxFlux-minFlux);
                if range<tol
                    range = 0;
                end
                FluxRange = [FluxRange; range];
            end   
        end
    end
    
    if ~isempty(FluxRange)
        FluxRange = sum(FluxRange);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = setBounds(model,indexes,forward,bounds)
    if length(indexes)>1
        if forward 
            %Sets an upper bound for backward reaction to avoid artificially 
            %induced high variability
            model.ub(indexes(2)) = bounds(2);
        else
            %Sets an upper bound for forward reaction with the maximum flux 
            %obtained in the previous step to avoid artificially induced
            %high variability
            model.ub(indexes(1)) = bounds(1);
        end
    end
end
