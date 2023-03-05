function [model, tunedKcats] = sensitivityTuning(model, desiredGrowthRate, modelAdapter, foldChange, verbose)
% sensitivityTuning
%    Function that relaxes the most limiting kcats until a certain growth rate
%    is reached. The function will update kcats in model.ec.kcat.
%
% Input:
%   model               an ecModel in GECKO 3 format (with ecModel.ec structure)
%   desiredGrowthRate  kcats will be relaxed until this growth rate is reached
%   modelAdapter       a loaded model adapter (Optional, will otherwise use the
%                      default model adapter).
%   foldChange         kcat values will be increased by this fold-change.
%                      (Opt, default 10)
%   verbose            logical whether progress should be reported (Optional,
%                      default true)
%
% Output:
%   model              ecModel with updated model.ec.kcat
%   tunedKcats         structure with information on tuned kcat values
%                      rxns     identifiers of reactions with tuned kcat
%                               values
%                      rxnNames names of the reactions in tunedKcats.rxns
%                      enzymes  enzymes that catalyze the reactions in
%                               tunedKcats.rxns, whose kcat value has been
%                               tuned.
%                      oldKcat  kcat values in the input model
%                      newKcat  kcat values in the output model, after tuning
%
% Usage:
%   [model, tunedKcats] = sensitivityTuning(model, desiredGrowthRate, modelAdapter, foldChange, verbose)

if nargin < 5 || isempty(verbose)
    verbose = true;
end
if nargin < 4 || isempty(foldChange)
    foldChange = 10;
end
if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
if nargin < 2 || isempty(desiredGrowthRate)
    desiredGrowthRate = params.gR_exp;
end

kcatList = [];
m = model;
m.c = double(strcmp(m.rxns, params.bioRxn));% Make sure that growth is maximized

[~,hs] = solveLP(m);
lastGrowth = 0;
if ~m.ec.geckoLight
    %for the full model, we first find the draw reaction with the most flux
    drawRxns = startsWith(m.rxns, 'usage_prot_');
    iteration = 1;
    while true
        [res,hs] = solveLP(m,0,[],hs); %skip parsimonius, only takes time
        if (lastGrowth == -res.f)
            disp('No growth increase from increased kcats - check if the constraints on the uptake reactions are too tight!')
            break;
        end
        lastGrowth = -res.f;
        if (lastGrowth >= desiredGrowthRate)
            break;
        end
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        iteration            = iteration + 1;
        %find the highest draw_prot rxn flux
        drawFluxes           = zeros(length(drawRxns),1);
        drawFluxes(drawRxns) = res.x(drawRxns);
        [~,sel]              = min(drawFluxes); % since bounds -1000 to 0
        %Now get the metabolite
        metSel               = m.S(:,sel) < 0; % negative coeff
        %now find the reaction with the largest consumption of this protein
        protFluxes           = m.S(metSel,:).' .* res.x; %negative
        [~,rxnSel]           = min(protFluxes);
        kcatList             = [kcatList, rxnSel];
        rxn                  = m.rxns(rxnSel);
        targetSubRxn         = strcmp(m.ec.rxns, rxn);
        m.ec.kcat(targetSubRxn) = m.ec.kcat(targetSubRxn) .* foldChange;
        m                    = applyKcatConstraints(m,targetSubRxn);
    end

else
    origRxns = extractAfter(m.ec.rxns,4);
    iteration = 1;
    while true
        res = solveLP(m,0); %skip parsimonius, only takes time
        if (lastGrowth == -res.f)
            disp('No growth increase from increased kcats - check if the constraints on the uptake reactions are too tight!')
            break;
        end
        lastGrowth = -res.f;
        if (lastGrowth >= desiredGrowthRate)
            break;
        end
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        iteration       = iteration + 1;
        %find the highest protein usage flux
        protPoolStoich  = m.S(strcmp(m.mets, 'prot_pool'),:).';
        [~,sel]         = min(res.x .* protPoolStoich); %max consumption
        kcatList        = [kcatList, sel];
        rxn             = m.rxns(sel.');
        targetSubRxns   = strcmp(origRxns, rxn);
        m.ec.kcat(targetSubRxns) = m.ec.kcat(targetSubRxns) .* foldChange;
        m               = applyKcatConstraints(m,rxn);
    end
end

kcatList        = unique(kcatList);
tunedKcats.rxns = m.rxns(kcatList);
tunedKcats.rxnNames = m.rxnNames(kcatList);
if ~model.ec.geckoLight
    [~, rxnIdx]     = ismember(tunedKcats.rxns,m.ec.rxns);
else
    [~, rxnIdx]     = ismember(tunedKcats.rxns,origRxns);
end
tunedKcats.enzymes  = cell(numel(kcatList),1);
for i=1:numel(rxnIdx)
    [~, metIdx]     = find(m.ec.rxnEnzMat(rxnIdx(i),:));
    tunedKcats.enzymes{i}   = strjoin(m.ec.enzymes(metIdx),';');
end
tunedKcats.oldKcat  = model.ec.kcat(rxnIdx);
tunedKcats.newKcat  = m.ec.kcat(rxnIdx);

model = m;
end
