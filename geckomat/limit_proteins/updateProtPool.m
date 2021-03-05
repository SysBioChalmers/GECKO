function model  = updateProtPool(model,Ptot,Pmeas,f,sigma)
%Calculate total mass of bounded enzymes
if any(isnan(model.concs))
    poolIndex  = find(strcmpi(model.rxns,'prot_pool_exchange'));
    difference = (Ptot*f-Pmeas);
    tempModel  = model;
    if (difference)>0
        %remove all measured enzymes mass (including flexibilizations) from
        %remaining pool
        newPool   = difference*sigma;
        tempModel = setParam(tempModel,'ub',poolIndex,newPool);
        solution  = solveLP(tempModel);
        if ~isempty(solution.x)
            fprintf('\necModel succesfully constrained with enzyme abundances\n\n');
            model = tempModel;
        else
            %Allow any saturation level
            sigma     = 1;
            newPool   = difference*sigma;
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
                fprintf('\necModel succesfully constrained with enzyme abundances\n\n')
                warning(['Sigma has been readjusted to: ' num2str(sigma)])
            else
                %If the optimal sigma factor is higher than one, then
                %calculate minimal protein pool requirements
                tempModel = setParam(tempModel,'ub',poolIndex,1000);
                solution  = solveLP(tempModel);
                PpoolUB   = solution.x(poolIndex);
                f         = (PpoolUB+Pmeas)/Ptot;
                fprintf(['\nf factor has been refitted to a value of ' num2str(f) ' according to\nminimal protein requirements of the constrained model.\n\n'])
                if f<=1
                    model = setParam(tempModel,'ub',poolIndex,1.01*PpoolUB);
                else
                    model = setParam(tempModel,'ub',poolIndex,1000);
                end
            end
            model.c = originalC;
        end
    elseif (difference)<0
        warning('The total measured protein mass exceeds the remaining protein content')
        warning('Returning an empty ecModel')
        model = [];
    else
         warning('The total measured protein mass is equal to the total protein content, therefore no protein pool is available for the unmeasured model enzymes')
         warning('Setting UB = 0 for prot_pool_exchange reaction')
         model = setParam(tempModel,'ub',poolIndex,0);
         sol   = solveLP(model);
         if isempty(sol.x)
             warning('The resulting ecModel is unfeasible under the imposed constraints')
         end
    end
end
end