function [model, proteins, frequence] = flexibilizeProtConcs(model, expGrowth, foldChange, iterationsPerEnzyme, modelAdapter)
% flexibilizeProtConcs
%   Flexibilize protein concentration of an ecModel with constrained with
%   proteomcis data
%
% Input:
%   model                 an ecModel in GECKO 3 version
%   expGrowth             Estimaed experimental growth rate. If not specified,
%                         the value will be read from the model adapter.
%   foldChange            a value how much increase the protein concentration.
%                         (Optional, default = 0.5)
%   iterationsPerEnzyme   the number of iterations that an enzyme can be increased.
%                         A zero number can be defined. if zero is defined no limit
%                         will be set, and it will increase the protein concentration
%                         until reach de defined growth rate (Optional, default = 5)
%   modelAdapter          a loaded model adapter (Optional, will otherwise use the
%                         default model adapter).
%
% Output:
%   model                 ecModel where the UB of measured protein have been increased
%                         to allow reach a defined growth rate.
%   proteins              a vector array with proteins evaluated
%   frequence             a vector array with a number of n times the UB was increase.
%                         Then new UB is [Ei] * 1.foldChange^n
%
% Usage:
%    [model, proteins, frequence] = flexibilizeProtConcs(model, expGrowth, foldChange, iterationsPerEnzyme, modelAdapter);

if nargin < 5 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if nargin < 4 || isempty(iterationsPerEnzyme)
    iterationsPerEnzyme = 5;
end

if nargin < 3 || isempty(foldChange)
    foldChange = 0.5;
end

if nargin < 2 || isempty(expGrowth)
    expGrowth = modelAdapter.getParameters().gR_exp;
end

% If a zero value is defined, not iteration limit will be set
if iterationsPerEnzyme == 0
    iterationsPerEnzyme = inf;
end

% In case the model have not been protein constrained
% model = constrainProtConcs(model);

sol = solveLP(model);
predGrowth = abs(sol.f);

% Get those proteins with a concentration defined
protConcs = find(~isnan(model.ec.concs));

% Get enzymes names with a concentration value
proteins = model.ec.enzymes(protConcs);

% Store how many time is predicted a enzyme with the highest coeff
frequence = zeros(numel(proteins),1);

if any(protConcs)
    while predGrowth < expGrowth

        % Get the control coefficients
        [~, controlCoeffs] = getConcControlCoeffs(model, proteins, foldChange, 0.75);

        % Find the maximum coefficient index
        [~,maxIdx] = max(controlCoeffs);

        % Increase +1 the frequence
        frequence(maxIdx) = frequence(maxIdx) + 1;

        % Allow to increase the UB maximum five times
        if frequence(maxIdx) <= iterationsPerEnzyme

            % Set how much will increase the UB
            increase = foldChange*frequence(maxIdx);

            % Get usage rxn for the highest coeff
            protUsageIdx = strcmpi(model.rxns, strcat('usage_prot_', proteins(maxIdx)));

            % replace the UB
            model.ub(protUsageIdx) = model.ec.concs(protConcs(maxIdx)) * (1+increase);

            % Get the new growth rate
            sol = solveLP(model);
            predGrowth = abs(sol.f);

            disp(['Protein ' proteins{maxIdx} ' UB adjusted. Grow: ' num2str(predGrowth)])
        else
            disp( ['Limit have been reached. Protein '  proteins{maxIdx} ' seems to be problematic. Consider changing the kcat '])
            break
        end
    end
else
    error('Protein concentrations have not been defined. Please run readProteomics and constrainProtConcs')
end

