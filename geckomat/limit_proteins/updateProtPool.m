function model  = updateProtPool(model,Ptot,Pmeas,f,sigma)
%Calculate total mass of bounded enzymes
if any(isnan(model.concs))
    poolIndex  = find(strcmpi(model.rxns,'prot_pool_exchange'));
    difference = Ptot-Pmeas;
    tempModel  = model;
    if (difference)>=0
        %remove all measured enzymes mass (including flexibilizations) from
        %remaining pool
        newPool   = (Ptot*f-Pmeas)*sigma;
        tempModel = setParam(tempModel,'ub',poolIndex,newPool);
        solution  = solveLP(tempModel);
        if ~isempty(solution.x)
            disp('ecModel succesfully constrained with enzyme abundances')
            model = tempModel;
        else
            %Allow any saturation level
            sigma     = 1;
            newPool   = (Ptot*f-Pmeas)*sigma;
            tempModel = setParam(tempModel,'ub',poolIndex,newPool);
            %As bounds for CUR and growth have already been set, protein
            %pool minimization is set as an objective
            originalC = model.c;
            tempModel = setParam(tempModel,'obj',poolIndex,-1);
            solution  = solveLP(tempModel);
            if ~isempty(solution.x)
                sigma = solution.x(poolIndex)/newPool;
                %Set minimal prot pool requirements as UB (plus a 0.1% of
                %flexibilization to avoid overconstraining)
                model = setParam(tempModel,'ub',poolIndex,1.001*solution.x(poolIndex));
                disp('ecModel succesfully constrained with enzyme abundances')
                warning(['Sigma has been readjusted to: ' num2str(sigma)])
            else
                %If the optimal sigma factor is higher than one, then protein
                %pool is unbounded
                warning('Unfeasible protein flexibilization, unconstraining protein pool for obtention of a feasible model')
                model = setParam(tempModel,'ub',poolIndex,1000);
            end
            model.c = originalC;
        end
    else
        warning('The total measured protein mass exceeds the total protein content.')
        model = [];
    end
end
end