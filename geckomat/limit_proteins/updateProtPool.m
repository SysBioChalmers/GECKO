function model  = updateProtPool(model,Ptot,fs)
%Calculate total mass of bounded enzymes
if any(isnan(model.concs))
    poolIndex  = find(strcmpi(model.rxns,'prot_pool_exchange'));
    Pmeasured  = sum(model.concs(~isnan(model.concs)));
    difference = Ptot-Pmeasured;
    tempModel  = model;
    if (difference)>=0
        %remove all measured enzymes mass (including flexibilizations) from
        %remaining pool
        newPool   = (difference)*fs;
        tempModel = setParam(tempModel,'ub',poolIndex,newPool);
        solution  = solveLP(tempModel);
        if ~isempty(solution.x)
            disp('ecModel succesfully constrained with enzyme abundances')
            model = tempModel;
        else
            %Allow any saturation level
            newPool   = difference;
            tempModel = setParam(tempModel,'ub',poolIndex,newPool);
            %As bounds for CUR and growth have already been set, protein
            %pool minimization is set as an objective
            originalC = model.c;
            tempModel = setParam(tempModel,'obj',poolIndex,-1);
            solution  = solveLP(tempModel);
            if ~isempty(solution.x)
                fs = solution.x(poolIndex)/newPool;
                model = setParam(tempModel,'ub',poolIndex,1.001*solution.x(poolIndex));
                disp('ecModel succesfully constrained with enzyme abundances')
                warning(['fs has been readjusted to: ' num2str(fs)])
            else
                %If the optimal fs factor is higher than one, then protein
                %pool is unbounded
                warning('Unfeasible protein flexibilization')
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