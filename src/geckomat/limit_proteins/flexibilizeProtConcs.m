function [model, flexProt] = flexibilizeProtConcs(model, expGrowth, foldChange, iterPerEnzyme, modelAdapter, verbose)
% flexibilizeProtConcs
%   Flexibilize protein concentration of an ecModel with constrained with
%   proteomics data. The upper bound of the protein usage reaction is
%   changed, while the concentrations in ecModel.ec.concs remain unchanged.
%
%   If no (more) limiting enzyme concentrations can be found, it might be
%   the protein pool exchange that is limiting growth. In that case, an
%   attempt will be made to relax the protein pool exchange reaction, and
%   if the growth rate indeed increases, it is suggested to set the protein
%   pool exchange unconstrained (lb=-1000) before again running
%   flexibilizeProtConcs. Such situations, where a proteomics integrated
%   ecModel is overconstrained, may occur if the ecModel should be able to
%   simulate maximum growth rate (from e.g. batch cultivation). 
%   
%   If relaxing the protein pool exchange does not increase the growth rate, then this
%   is not due to enzyme constraints, but rather an issue with the
%   metabolic network itself, or the set nutrient exchange is not
%   sufficient.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   expGrowth       estimated experimental growth rate. If not specified,
%                   the value will be read from the model adapter
%   foldChange      a value how much increase the protein concentration
%                   (Optional, default = 2)
%   iterPerEnzyme   the number of iterations that an enzyme can be increased.
%                   A zero number can be defined. if zero is defined no limit
%                   will be set, and it will increase the protein concentration
%                   until reach de defined growth rate (Optional, default = 5)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter)
%   verbose         logical whether progress should be reported (Optional,
%                   default true)
%
% Output:
%   model           ecModel where the UB of measured protein have been increased
%                   to allow reach a defined growth rate
%   flexProt        array with information about flexibilized proteins
%                   uniprotIDs  enzymes whose usage UB was flexibilized
%                   oldConcs    original concentrations, from mode.ec.concs
%                   flexConcs   flexibilized concentrations, new UB in
%                               model
%                   ratioIncr   ratio by which the concentration increased,
%                               the enzymes will be sorted by this field
%                   frequence   numeric how often the enzyme has been
%                               step-wise flexibilized
%
% Usage:
%   [model, flexProt] = flexibilizeProtConcs(model, expGrowth, foldChange, iterPerEnzyme, modelAdapter, verbose)

if nargin < 6 || isempty(verbose)
    verbose = true;
end
if nargin < 5 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
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

bioRxnIdx = getIndexes(model,params.bioRxn,'rxns');
if model.ub(bioRxnIdx) < expGrowth
    fprintf(['The upper bound of the biomass reaction was %s, while expGrowth '...
             'is %s. The upper bound has been set to expGrowth.'],...
             num2str(model.ub(bioRxnIdx)), num2str(expGrowth))
    model.ub(bioRxnIdx) = expGrowth;
end

% Keep track if protein pool will be flexibilized
changedProtPool = false;

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

        % Stop looking for limiting proteins all coeffs are zero
        if any(controlCoeffs)

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
                model.lb(protUsageIdx) = -(model.ec.concs(protConcs(maxIdx)) * (1+increase));

                % Get the new growth rate
                sol = solveLP(model,0,[],hs);
                predGrowth = abs(sol.f);

                if verbose
                    disp(['Protein ' proteins{maxIdx} ' LB adjusted. Grow: ' num2str(predGrowth)])
                end
            else
                disp(['Limit has been reached. Protein '  proteins{maxIdx} ' seems to be problematic. Consider changing the kcat.'])
                flexBreak=true;
                break
            end
        else
            disp('No (more) limiting proteins have been found. Attempt to increase protein pool exchange.')
            protPoolIdx = getIndexes(model,'prot_pool_exchange','rxns');
            oldProtPool = -model.lb(protPoolIdx);
            tempModel = setParam(model,'lb',protPoolIdx,-1000);
            tempModel = setParam(tempModel,'ub',bioRxnIdx,expGrowth);
            sol = solveLP(tempModel);
            if (abs(sol.f)-predGrowth)>1e-10 % There is improvement in growth rate
                % Find new protein pool constraint
                predGrowth = abs(sol.f);
                tempModel = setParam(tempModel,'lb',bioRxnIdx,predGrowth);
                tempModel = setParam(tempModel,'obj',protPoolIdx,1);
                sol = solveLP(tempModel);
                newProtPool = abs(sol.f);
                model.lb(protPoolIdx) = -newProtPool;
                fprintf(['Changing the lower bound of protein pool exchange from %s ',...
                         'to %s enabled a\ngrowth rate of %s. It can be helpful to set ', ...
                         'the lower bound of protein\npool exchange to -1000 before ',...
                         'running flexibilizeProtConcs.\n'],num2str(-oldProtPool), ...
                         num2str(-newProtPool),num2str(predGrowth));
                changedProtPool = true;
            else
                fprintf(['Protein pool exchange not limiting. Inability to reach growth ',...
                        'rate is not related to\nenzyme constraints. Maximum growth rate ',...
                        'is %s.\n'], num2str(predGrowth))
            end
            % To return an usable model with these conditions, set the
            % highest growth rate reached as the experimental value
            expGrowth = predGrowth;
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

if flexBreak == false && ~isempty(protFlex)
    % Not all flexibilized proteins require flexibilization in the end
    % Test flux distribution with minimum prot_pool usage
    modelTemp       = setParam(model,'var',params.bioRxn,expGrowth,0.5);
    modelTemp       = setParam(modelTemp,'obj','prot_pool_exchange',1);
    sol             = solveLP(modelTemp,1);
    % Check which concentrations can remain unchanged
    newConcs        = abs(sol.x(protUsageIdx));
    keepOld         = newConcs<oldConcs;
    % Set UBs
    model.lb(protUsageIdx(keepOld))  = -oldConcs(keepOld);
    model.lb(protUsageIdx(~keepOld)) = -newConcs(~keepOld);
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
    flexProt.flexConcs  = abs(model.lb(protUsageIdx));
    flexProt.frequence  = frequence(frequence>0);
end
if changedProtPool
    flexProt.uniprotIDs(end+1) = {'prot_pool'};
    flexProt.oldConcs(end+1)   = oldProtPool;
    flexProt.flexConcs(end+1)  = newProtPool;
    flexProt.frequence(end+1)  = 1;
end
% Sort by biggest ratio increase
flexProt.ratioIncr  = flexProt.flexConcs./flexProt.oldConcs;
[~,b] = sort(flexProt.ratioIncr,'descend');
flexProt.ratioIncr  = flexProt.ratioIncr(b);
flexProt.uniprotIDs = flexProt.uniprotIDs(b);
flexProt.oldConcs   = flexProt.oldConcs(b);
flexProt.flexConcs  = flexProt.flexConcs(b);
flexProt.frequence  = flexProt.frequence(b);
end
