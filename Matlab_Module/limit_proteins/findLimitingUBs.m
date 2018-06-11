function [maxIndex,flag] = findLimitingUBs(model,protIndxs,flexFactor,option)
% findLimitingUBs
%   Find the top growth-limiting enzyme or enzyme complex in a given ecModel
%   overconstrained with proteomics data.
%   
%   model       ecModel with proteomic constraints (individual enzyme
%               levels)
%   protIndxs   Indexes of the protein exchange reactions (measured
%               enzymes)
%   flexFactor  Flexibilization factor used to relax the upper bounds in
%               every iteration
%   option      1 if the sensitivity analysis is done on the individual
%               enzyme exchanges or 2, if it is done on a reaction basis, 
%               for cases in which the limitation might be an enzyme c
%               complex rather than an individual enzyme.
%
%   maxIndex    Index of the corresponding limiting enzyme exchange
%               reaction, this variable has a length>1 for cases in which 
%               the limitation was found to be a complex, every index 
%               correspond to the exchange reaction for every subunit in 
%               the complex.
%   flag        TRUE if a limiting enzyme or complex was found
%
%   Usage: [maxIndex,flag] = findLimitingUBs(model,protIndxs,flexFactor,option)
%
%   Ivan Domenzain, 2018-06-11
%

%Find objective reaction and perform an initial simulation for getting the
%base flux distribution
objIndex  = find(model.c==1);
solution  = solveLP(model,1);
tolerance = 1E-2;
flag      = false;

switch option
    %Analyse enzyme usages upper bounds
    case 1
        coeffs = zeros(length(protIndxs),1);
        % flexibilize ub for every protein exchange in a temporal model
        for i=1:length(protIndxs)
            index = protIndxs(i);
            %The analyis is run just on those enzymes which usage is equal or
            %very close to its respective upper bound 
            if ((model.ub(index)-solution.x(index))/model.ub(index))<=tolerance
                tempModel           = model;
                tempModel.ub(index) = tempModel.ub(index)*flexFactor;
                %get a new solution with the flexibilized ub and calculate the
                %effect on the growth rate prediction
                newSol = solveLP(tempModel,1);
                if ~isempty(newSol.f)
                    deltaUsage = newSol.x(index)-solution.x(index);
                    if deltaUsage~=0
                        deltaGrowth  = newSol.x(objIndex)-solution.x(objIndex);
                        ControlCoeff = deltaGrowth/deltaUsage;
                        coeffs(i)    = abs(ControlCoeff);
                    end
                end
            end
        end
        [~,maxIndex] = max(coeffs);
        if coeffs(maxIndex)~=0
            flag = true;
        end
        maxIndex = protIndxs(maxIndex);
    
    %Analyse reaction-associated enzymes UBs (for enzyme complexes)
    case 2
       solution = solution.x;
       coeffs   = zeros(length(model.rxns),1); 
       protRxns = find(contains(model.rxnNames,'prot_'));
       for i=1:protRxns(1)%length(model.rxns)
           %for flux carrying reactions
           if solution(i)~=0
               tempModel   = model;
               measurement = false;
               protIndexes = extractProteinExchanges(i,model);
               % If there are measured proteins associated with the i-th
               % rxn their respective upper bounds are flexibilized 
               for j = 1:length(protIndexes)
                   pos     = cell2mat(protIndexes(j));
                   relDiff = (model.ub(pos)-solution(pos))/model.ub(pos);
                   if model.ub(pos) < 1000 && relDiff <= tolerance
                       tempModel.ub(pos) = tempModel.ub(pos)*flexFactor;
                       measurement       = true;
                   end
               end
               %If proteins were flexibilized for the i-th reaction, a new
               %solution is obtained and the objective function sensitivity
               %is calculated
               if measurement
                   new_sol  = solveLP(tempModel,1);
                   deltaRxn = new_sol.x(i)-solution(i);
                    if deltaRxn~=0
                        deltaGrowth  = new_sol.x(objIndex)-solution(objIndex);
                        ControlCoeff = deltaGrowth/deltaRxn;
                        coeffs(i)    = abs(ControlCoeff);
                    end
               end
           end
       end
       [~,maxIndex] = max(coeffs);
        if coeffs(maxIndex)~=0
            flag = true;
        end
        maxIndex = extractProteinExchanges(maxIndex,model);
        maxIndex = cell2mat(maxIndex);
end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the indexes of protein exchange reactions associated to a metabolic
%reaction
function indexes = extractProteinExchanges(rxn,model)
indexes   = {};
rxn_mets  = model.mets(full(model.S(:,rxn)) < 0);
rxn_prots = find(~cellfun(@isempty,strfind(rxn_mets,'prot_')));
% If there are measured proteins associated with the rxn their respective
% upper bounds are flexibilized
for j = 1:length(rxn_prots)
    pos     = find(strcmpi([rxn_mets{rxn_prots(j)} '_exchange'],model.rxns));
    if ~isempty(pos)
        indexes = [indexes; pos];
    end
end
end