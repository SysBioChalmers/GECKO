function [model, flexProt] = flexibilizeProtConcs(model, expGrowth, foldChange, iterPerEnzyme, modelAdapter, verbose)
% flexibilizeProtConcs
%   Flexibilize protein concentration of an ecModel with constrained with
%   proteomics data. The upper bound of the protein usage reaction is
%   changed, while the 
%
% Input:
%   model           an ecModel in GECKO 3 version
%   expGrowth       Estimaed experimental growth rate. If not specified,
%                   the value will be read from the model adapter.
%   foldChange      a value how much increase the protein concentration.
%                   (Optional, default = 2)
%   iterPerEnzyme   the number of iterations that an enzyme can be increased.
%                   A zero number can be defined. if zero is defined no limit
%                   will be set, and it will increase the protein concentration
%                   until reach de defined growth rate (Optional, default = 5)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%   verbose         logical whether progress should be reported (Optional,
%                   default true)
%
% Output:
%   model           ecModel where the UB of measured protein have been increased
%                   to allow reach a defined growth rate.
%   flexProt        array with information about flexibilized proteins
%                   uniprotIDs  enzymes whose usage UB was flexibilized
%                   oldConcs    original concentrations, from mode.ec.concs
%                   flexConcs   flexibilized concentrations, new UB in
%                               model
%                   frequence   numeric how often the enzyme has been
%                               step-wise flexibilized
%
% Usage:
%    [model, flexProt] = flexibilizeProtConcs(model, expGrowth, foldChange, iterPerEnzyme, modelAdapter)
if nargin < 6 || isempty(verbose)
    verbose = true;
end
if nargin < 5 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 4 || isempty(iterPerEnzyme)
    iterPerEnzyme = 5;
end

if nargin < 3 || isempty(foldChange)
    foldChange = 2;
end

if nargin < 2 || isempty(expGrowth)
    expGrowth = modelAdapter.getParameters().gR_exp;
end

% If a zero value is defined, not iteration limit will be set
if iterPerEnzyme == 0
    iterPerEnzyme = inf;
end

% In case the model have not been protein constrained
% model = constrainProtConcs(model);

[sol,hs] = solveLP(model);
predGrowth = abs(sol.f);

% Get those proteins with a concentration defined
protConcs = find(~isnan(model.ec.concs));

% Get enzymes names with a concentration value
proteins = model.ec.enzymes(protConcs);

% Store how many time is predicted a enzyme with the highest coeff
frequence = zeros(numel(proteins),1);

flexBreak=false;
if any(protConcs)
    while predGrowth < expGrowth

        % Get the control coefficients
        [~, controlCoeffs] = getConcControlCoeffs(model, proteins, foldChange, 0.75);

        % Find the maximum coefficient index
        [~,maxIdx] = max(controlCoeffs);

        % Increase +1 the frequence
        frequence(maxIdx) = frequence(maxIdx) + 1;

        % Allow to increase the UB maximum five times
        if frequence(maxIdx) <= iterPerEnzyme

            % Set how much will increase the UB
            increase = foldChange*frequence(maxIdx);

            % Get usage rxn for the highest coeff
            protUsageIdx = strcmpi(model.rxns, strcat('usage_prot_', proteins(maxIdx)));

            % replace the UB
            model.ub(protUsageIdx) = model.ec.concs(protConcs(maxIdx)) * (1+increase);

            % Get the new growth rate
            sol = solveLP(model,0,[],hs);
            predGrowth = abs(sol.f);

            if verbose
                disp(['Protein ' proteins{maxIdx} ' UB adjusted. Grow: ' num2str(predGrowth)])
            end
        else
            disp(['Limit has been reached. Protein '  proteins{maxIdx} ' seems to be problematic. Consider changing the kcat '])
            flexBreak=true;
            break
        end
    end
else
    error('Protein concentrations have not been defined. Please run readProteomics and constrainProtConcs')
end

protFlex        = proteins(frequence>0);
protUsageIdx    = find(ismember(model.rxns, strcat('usage_prot_', protFlex)));
ecProtId        = find(ismember(model.ec.enzymes,protFlex));
oldConcs        = model.ec.concs(ecProtId);

if flexBreak == false
    % Not all flexibilized proteins require flexibilization in the end
    % Test flux distribution with minimum prot_pool usage
    modelTemp       = setParam(model,'eq',params.bioRxn,expGrowth);
    modelTemp       = setParam(modelTemp,'obj','prot_pool_exchange',-1);
    sol             = solveLP(modelTemp,1,[],hs);
    % Check which concentrations can remain unchanged
    newConcs        = sol.x(protUsageIdx);
    keepOld         = newConcs<oldConcs;
    % Set UBs
    model.ub(protUsageIdx(keepOld))  = oldConcs(keepOld);
    model.ub(protUsageIdx(~keepOld)) = newConcs(~keepOld);
    % Output structure
    protFlex(keepOld) = [];
    oldConcs(keepOld) = [];
    newConcs(keepOld) = [];
    frequence=frequence(frequence>0);
    frequence(keepOld) = [];
    
    flexProt.uniprotIDs = protFlex;
    flexProt.oldConcs   = oldConcs;
    flexProt.flexConcs  = newConcs;
    flexProt.frequence  = frequence(frequence>0);
else
    protFlex        = proteins(frequence>0);
    protUsageIdx    = find(ismember(model.rxns, strcat('usage_prot_', protFlex)));

    flexProt.uniprotIDs = protFlex;
    flexProt.oldConcs   = model.ec.concs(ecProtId);
    flexProt.flexConcs  = model.ub(protUsageIdx);
    flexProt.frequence  = frequence(frequence>0);
end
