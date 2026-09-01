function [enz, controlCoeffs] = getConcControlCoeffs(model, varargin)
% getConcControlCoeffs  Calculate control coefficients of protein usage.
%
% For each candidate protein whose usage reaction is running above the
% given limit fraction of its upper bound, temporarily increases that
% protein's usage upper bound by foldChange, re-solves, and expresses the
% resulting growth-rate increase per unit of concentration increase as a
% control coefficient. Proteins below the limit, or whose increase does
% not measurably raise growth, get a coefficient of zero.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% proteins : cell
%     a list of proteins to calculate the coefficients (default
%     model.ec.enzymes).
% foldChange : double
%     a value how much to increase the protein concentration (default 2).
% limit : double
%     threshold on the usage/concentration ratio (how saturated a
%     protein's usage upper bound already is); only proteins with a ratio
%     above this value are evaluated for a control coefficient (default
%     0, i.e. any nonzero usage is evaluated).
%
% Returns
% -------
% enz : logical
%     flags which of the input proteins had a usage/concentration ratio
%     above limit and were evaluated for a control coefficient.
% controlCoeffs : double
%     control coefficient per input protein (0 for proteins not evaluated
%     or whose increase did not measurably raise growth).
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [enz, controlCoeffs] = getConcControlCoeffs(model);
%     [enz, controlCoeffs] = getConcControlCoeffs(model, 'foldChange', 3);
%
% See also
% --------
% flexibilizeEnzConcs

p = parseGECKOargs(varargin, { ...
    'proteins',   []; ...
    'foldChange', []; ...
    'limit',      []});
proteins   = p.proteins;
foldChange = p.foldChange;
limit      = p.limit;

if isempty(limit)
    limit = 0;
end

if isempty(foldChange)
    foldChange = 2;
end

if isempty(proteins)
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
    prevConc = model.ub(protUsageRxnIdx(i));

    % Only consider those with a usage close the UB
    if (sol.x(protUsageRxnIdx(i))/prevConc) > limit
        
        % Update the logical vector
        enz(i) = 1;

        % Create a temporal model since coeff will be calculated one enzyme at
        % the time, without other change
        tempModel = model;
        % Increase the concentration by flexfactor
        newConc = prevConc*(foldChange);
        tempModel.ub(protUsageRxnIdx(i)) = newConc;

        % Get the new growth rate after the adjustment
        [tempSol,hs] = solveLP(tempModel,0,[],hs);
        tempGrowth = tempSol.f;
        
        % Calculate the coeff only if new growth rate is significantly
        % higher than initial value
        if (tempGrowth-initialGrowth)>1e-10
            controlCoeffs(i) = (tempGrowth-initialGrowth)/(newConc-prevConc);
        end
    end

end
end

