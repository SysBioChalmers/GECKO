function [model, tunedKcats] = sensitivityTuning(model, varargin)
% sensitivityTuning  Relax the most limiting kcats to reach a growth rate.
%
% Relaxes the most limiting kcats until a certain growth rate is reached. The
% function will update kcats in model.ec.kcat.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% desiredGrowthRate : double
%     kcats will be relaxed until this growth rate is reached (default: the
%     experimental growth rate params.gR_exp from the model adapter).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% foldChange : double
%     kcat values will be increased by this fold-change (default 10).
% protToIgnore : cell
%     vector of protein ids to be ignored in tuned kcats, e.g. {'P38122',
%     'Q99271'} (default []).
% verbose : logical
%     whether progress should be reported (default true).
%
% Returns
% -------
% model : struct
%     ecModel with updated model.ec.kcat.
% tunedKcats : struct
%     structure with information on tuned kcat values.
%
% Notes
% -----
% The tunedKcats structure has the following fields:
%
% - rxns : identifiers of reactions with tuned kcat values.
% - rxnNames : names of the reactions in tunedKcats.rxns.
% - enzymes : enzymes that catalyze the reactions in tunedKcats.rxns, whose
%   kcat value has been tuned.
% - oldKcat : kcat values in the input model.
% - newKcat : kcat values in the output model, after tuning.
%
% The model.ec.notes field will contain the original kcat value and source,
% unless the kcat has previously been set by sensitivityTuning, in which case
% the notes field remains unchanged.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, tunedKcats] = sensitivityTuning(model);
%     [model, tunedKcats] = sensitivityTuning(model, 'foldChange', 10);
%     [model, tunedKcats] = sensitivityTuning(model, desiredGrowthRate, modelAdapter, foldChange, protToIgnore, verbose);

p = parseGECKOargs(varargin, { ...
    'desiredGrowthRate', []; ...
    'modelAdapter',      []; ...
    'foldChange',        10; ...
    'protToIgnore',      {}; ...
    'verbose',           true});
desiredGrowthRate = p.desiredGrowthRate;
modelAdapter      = p.modelAdapter;
foldChange        = p.foldChange;
protToIgnore      = p.protToIgnore;
verbose           = p.verbose;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
if isempty(desiredGrowthRate)
    desiredGrowthRate = params.gR_exp;
end

kcatList = [];
m = model;
m.c = double(strcmp(m.rxns, params.bioRxn));% Make sure that growth is maximized

[res,hs] = solveLP(m);
if isempty(res.x)
    error('FBA of input model gives no valid result. Reduce protein pool constraint with setProtPoolSize and/or check if exchange constraints are correctly defined.')
end
lastGrowth = 0;
if ~m.ec.geckoLight
    %for the full model, we first find the draw reaction with the most flux
    drawRxns = startsWith(m.rxns, 'usage_prot_');
    idxToIgnore = cellfun(@(x) find(strcmpi(model.rxns, ['usage_prot_' x])), protToIgnore);
    iteration = 1;
    while true
        [res,hs] = solveLP(m,0,[],hs); %skip parsimonius, only takes time
        if (lastGrowth == res.f)
            printOrange('WARNING: No growth increase from increased kcats - check if the constraints on the uptake reactions are too tight.\n');
            break;
        end
        lastGrowth = res.f;
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        if (lastGrowth >= desiredGrowthRate)
            break;
        end
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        iteration            = iteration + 1;
        %find the highest draw_prot rxn flux
        drawFluxes           = zeros(length(drawRxns),1);
        drawFluxes(drawRxns) = res.x(drawRxns);
        % Remove from the list user defined proteins
        drawFluxes(idxToIgnore) = 0;
        [~,sel]              = max(drawFluxes); % since bounds 0 to 1000
        %Now get the metabolite
        metSel               = m.S(:,sel) > 0; % positive coeff (forward usage)
        %now find the reaction with the largest consumption of this protein
        protFluxes           = m.S(metSel,:).' .* res.x; %negative
        [~,rxnSel]           = min(protFluxes);
        kcatList             = [kcatList, rxnSel];
        rxn                  = m.rxns(rxnSel);
        targetSubRxn         = strcmp(m.ec.rxns, rxn);
        if ~strcmp(m.ec.source(targetSubRxn),'sensitivityTuning')
            oldNote          = m.ec.notes{targetSubRxn};
            newNote          = ['preTuneKcat=' num2str(m.ec.kcat(targetSubRxn)) ' | source:' m.ec.source{targetSubRxn}];
            if ~isempty(oldNote)
                newNote      = [oldNote '; ' newNote];
            end
            m.ec.notes{targetSubRxn}    = newNote;
        end
        m.ec.kcat(targetSubRxn)     = m.ec.kcat(targetSubRxn) .* foldChange;
        m.ec.source(targetSubRxn)   = {'sensitivityTuning'};
        m                    = applyKcatConstraints(m,targetSubRxn);
    end

else
    origRxns = extractAfter(m.ec.rxns,4);
    %find the reactions involved in proteins to be ignored
    idxToIgnore = cellfun(@(x) find(m.ec.rxnEnzMat(:, strcmpi(m.ec.enzymes, x))), protToIgnore, 'UniformOutput', false);
    %create an unique vector
    idxToIgnore = unique(cat(1, idxToIgnore{:}));
    %get the correct idx in model.rxns
    idxToIgnore = cellfun(@(x) find(strcmpi(m.rxns, x)), origRxns(idxToIgnore));
    iteration = 1;
    while true
        res = solveLP(m,0); %skip parsimonius, only takes time
        if (lastGrowth == res.f)
            printOrange('No growth increase from increased kcats - check if the constraints on the uptake reactions are too tight.\n');
            break;
        end
        lastGrowth = res.f;
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
        protPoolStoich(idxToIgnore) = 0;
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
tunedKcats.source   = model.ec.source(rxnIdx);

model = m;
end
