function FluxRange = MAXmin_Optimizer(model,indexes,FixedValues,tol)
% MAXmin_Optimizer
%
%   Get a model and the index(es) of the rxns to maximize and minimize. If
%   both optimizations were feasible, then a FV range is returned.
%
%   Usage: FluxRange = MAXmin_Optimizer(model,indexes,FixedValues,tol)
%
% Ivan Domenzain.      Last edited: 2019-12-04

FluxRange = [];
%Index is a 2 cells array when the rxn is a splitted reversible rxn
for i=1:length(indexes)
    %%% Maximization
    %Set objective function to maximization of the desired flux
    Temp_model = setParam(model,'obj',indexes(i),1);
    if length(indexes)>1
        if i==1
            %For optimization of forward rxn then backward rxn should
            %be blocked
            Temp_model.ub(indexes(2)) = FixedValues(2);
        else
            %For optimization of backward rxn then forward rxn should
            %be blocked
            Temp_model.ub(indexes(1)) = FixedValues(1);
        end
    end
    %Block the flux for the opposite rxn for irrev models
    solution = solveLP(Temp_model);
    %If Maximization was feasible, then proceed to minimization %%%%
    if ~isempty(solution.f)
        maxFlux    = solution.x(indexes(i));
        Temp_model = setParam(Temp_model,'obj',indexes(i),-1);
        solution   = solveLP(Temp_model);
        if ~isempty(solution.f)
            minFlux   = solution.x(indexes(i));
            FluxRange = [FluxRange; maxFlux-minFlux];
        end
    end
end
%Both rxns (irrev case) variability should be considered when
%calculating the overall rxn flux range
if ~isempty(FluxRange)
    FluxRange = sum(FluxRange);
    if FluxRange<tol
        FluxRange = 0;
    end
end
end