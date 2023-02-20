function model = flexibilizeProtConcs(model, foldChange, expGrowth, modelAdapter)
% flexibilizeProtConcs
%   Flexibilize protein concentration of an ecModel with constrained with
%   proteomcis data
%
% Input:
%   model           an ecModel in GECKO 3 version
%   foldChange      a value how much increase the protein concentration.
%                   (Optional, default = 0.1)
%   expGrowth       Estimaed experimental growth rate. If not specified,
%                   the value will be read from the model adapter.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
% Output:
%   enz             a logical vector of enzymes analyzed
%   controlCoeffs   a vector array with the coefficients
%
% Usage:
%    [enz, controlCoeffs] = getConcControlCoeffs(model, proteins, foldChange, limit);

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if nargin < 3 || isempty(expGrowth)
    expGrowth = modelAdapter.getParameters().gR_exp;
end

if nargin < 2
    foldChange = 0.1;
end

% In case the model have not been protein constrained
% model = constrainProtConcs(model);

sol = solveLP(model);
predGrowth = abs(sol.f);

% Get those proteins with a concentration defined
protConcs = ~isnan(model.ec.concs);

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
        if frequence(maxIdx) <= 5

            % Set how much will increase the UB
            increase = foldChange*frequence(maxIdx);

            % Get usage rxn for the highest coeff
            protUsageIdx = strcmpi(model.rxns, strcat('usage_prot_', proteins(maxIdx)));

            % replace the UB
            model.ub(protUsageIdx) = model.ub(protUsageIdx) * (1+increase);

            % Get the new growth rate
            sol = solveLP(model);
            predGrowth = abs(sol.f);

            disp(['Protein ' proteins{maxIdx} ' UB adjusted. Grow: ' num2str(predGrowth)])
        else
            disp( ['Protein '  proteins{maxIdx} ' seems to be problematic. Consider changing the kcat '])
            break
        end

    end
else
    error('Protein concentrations have not been defined. Please run readProteomics and constrainProtConcs')
end

