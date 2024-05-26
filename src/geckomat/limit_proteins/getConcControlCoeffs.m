function [enz, controlCoeffs] = getConcControlCoeffs(model, proteins, foldChange, limit)
% getConcControlCoeffs
%   Calculate an control coefficients of protein usage.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   proteins        a list of proteins to calculate the coefficients. (Optional,
%                   default = model.ec.enzymes)
%   foldChange      a value how much increase the protein concentration.
%                   (Optional, default = 2)
%   limit           a value to determine limiting protein usage reactions.
%                   Calculate as usage/concentration (Optional, default = 0)
%
% Output:
%   enz             a logical vector of enzymes analyzed
%   controlCoeffs   a vector array with the coefficients
%
% Usage:
%    [enz, controlCoeffs] = getConcControlCoeffs(model, proteins, foldChange, limit);

if nargin < 4 || isempty(limit)
    limit = 0;
end

if nargin < 3 || isempty(foldChange)
    foldChange = 2;
end

if nargin < 2 || isempty(proteins)
    proteins = model.ec.enzymes;
end

% for now 
enz = false(length(proteins),1);
controlCoeffs = zeros(length(proteins),1);

[sol,hs] = solveLP(model);
initialGrowth = sol.f;

% Get enzyme index
[~, protIdx] = ismember(proteins, model.ec.enzymes);

% Get the protein usage reactions
protUsageRxns = strcat('usage_prot_', model.ec.enzymes(protIdx));
[~, protUsageRxnIdx] = ismember(protUsageRxns, model.rxns);

for i = 1:numel(proteins)
    % Get the previous concentration
    prevConc = model.lb(protUsageRxnIdx(i));

    % Only consider those with a usage close the UB
    if (sol.x(protUsageRxnIdx(i))/prevConc) > limit
        
        % Update the logical vector
        enz(i) = 1;

        % Create a temporal model since coeff will be calculated one enzyme at
        % the time, without other change
        tempModel = model;
        % Increase the concentration by flexfactor
        newConc = prevConc*(foldChange);
        tempModel.lb(protUsageRxnIdx(i)) = newConc;

        % Get the new growth rate after the adjustment
        [tempSol,hs] = solveLP(tempModel,0,[],hs);
        tempGrowth = tempSol.f;
        
        % Calculate the coeff only if new growth rate is significantly
        % higher than initial value
        if (tempGrowth-initialGrowth)>1e-10
            controlCoeffs(i) = (tempGrowth-initialGrowth)/(prevConc-newConc);
        end
    end

end
end

