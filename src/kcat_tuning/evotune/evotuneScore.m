function [objective,rmse,rmseDetail,ecModel] = evotuneScore(ecModel,evoData,varargin)
% evotuneScore  Score one ecModel's current kcats against evotune data.
%
% Internal helper used by cmaesKcatTuning, screenKcatLeverage and
% tunePriorPenaltyWeight to score one ecModel (with a given set of kcats
% already applied via applyKcatConstraints) against the experimental
% growth and flux data. Fuses what geckopy splits into simulate.py (the
% FBA half) and distance.py (the RMSE half) into one function, matching
% this module's own predecessor abc_max.m -- there is no MATLAB test
% suite here to benefit from Python's testability-driven split.
%
% Parameters
% ----------
% ecModel : struct
%     ecModel with the current candidate kcats already applied (via
%     applyKcatConstraints). Not mutated.
% evoData : struct
%     structure with experimental data to be used, as returned by
%     loadEvotuneData.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% maxGrowthWeight : double
%     relative weight of the max-growth dataset against the flux dataset:
%     (rmseFlux + w*rmseMaxGrate) / (w+1) (default 1.0).
% priorPenaltyWeight : double
%     weight on a prior term added to the objective (default 0, i.e. the
%     objective equals the plain RMSE). See kcatVec/kcat0/sigma0log below.
% kcatVec, kcat0, sigma0log : double
%     only used when priorPenaltyWeight > 0: the tunable kcat subset's
%     current values, prior values, and per-kcat prior log-space std dev,
%     all the same length and in the same order. The penalty added to the
%     objective is priorPenaltyWeight * mean((log(kcatVec./kcat0)./sigma0log).^2).
%
% Returns
% -------
% objective : double
%     rmse, plus the prior penalty term when priorPenaltyWeight > 0. This
%     is what a search should minimise.
% rmse : double
%     the plain (unpenalised) average RMSE, always comparable across runs
%     regardless of priorPenaltyWeight.
% rmseDetail : cell
%     per-condition RMSE values, labelled by data type and condition index.
% ecModel : struct
%     the input ecModel, with ec.excarbon cached on it if it wasn't
%     already present. MATLAB structs are copied on assignment, not
%     shared by reference like a Python object, so callers scoring many
%     candidates in a loop should carry this output back in as the base
%     copy for the next candidate, rather than recomputing excarbon (an
%     expensive metabolite-formula scan) on every single call.
%
% See also
% --------
% cmaesKcatTuning, screenKcatLeverage, loadEvotuneData

p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'maxGrowthWeight', 1.0; ...
    'priorPenaltyWeight', 0.0; ...
    'kcatVec', []; ...
    'kcat0', []; ...
    'sigma0log', []});
modelAdapter        = p.modelAdapter;
maxGrowthWeight     = p.maxGrowthWeight;
priorPenaltyWeight  = p.priorPenaltyWeight;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

if ~isfield(ecModel,'excarbon')
    ecModel = computeExcarbon(ecModel);
    ecModel.excarbon(ecModel.excarbon == 0) = 1;
end

rmse_1 = []; rmseDetail1 = []; rmse_2 = []; rmseDetail2 = [];

if ~isempty(evoData.fluxData)
    [rmse_1, rmseDetail1] = simulateCondition(ecModel,evoData.fluxData,true,evoData.zeroFlux,modelAdapter);
    rmseDetail1 = [cellstr("fluxData_" + string(1:numel(rmseDetail1)))',num2cell(rmseDetail1)];
end
if ~isempty(evoData.maxGrate)
    [rmse_2, rmseDetail2] = simulateCondition(ecModel,evoData.maxGrate,false,evoData.zeroFlux,modelAdapter);
    rmseDetail2 = [cellstr("maxGrowth_" + string(1:numel(rmseDetail2)))',num2cell(rmseDetail2)];
end

parts = []; weights = [];
if ~isempty(rmse_1)
    parts(end+1) = rmse_1; weights(end+1) = 1.0;
end
if ~isempty(rmse_2)
    parts(end+1) = rmse_2; weights(end+1) = maxGrowthWeight;
end
if isempty(parts)
    rmse = NaN;
else
    rmse = sum(parts .* weights) / sum(weights);
end
rmseDetail = [rmseDetail1;rmseDetail2];

objective = rmse;
if priorPenaltyWeight > 0 && ~isempty(p.kcatVec) && ~isempty(p.kcat0) && ~isempty(p.sigma0log)
    dev = log(p.kcatVec(:) ./ p.kcat0(:)) ./ p.sigma0log(:);
    objective = rmse + priorPenaltyWeight * mean(dev.^2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rmse, rmseList] = simulateCondition(ecModel,data,constrain,rxn2block,modelAdapter)
% One condition's FBA simulation + carbon-weighted RMSE against measured
% data. Ported from abc_max.m's rmsecal.

rmseList = zeros(length(data.conds),1);

[~,exchIdx] = ismember(data.conds,data.exchMets);
if any(exchIdx == 0)
    missingMets = strjoin(unique(data.conds(exchIdx == 0)), '; ');
    error('Carbon source(s) "%s" in the provided fluxData or maxGrowth cannot be matched by name with an exchange reaction.', missingMets);
end

for i = 1:length(data.conds)
    % Set all other carbon sources to zero, then unblock this condition's
    model_tmp = setParam(ecModel,'lb',data.exchRxnIDs(unique(exchIdx)),0);
    if constrain % RMSE from flux data: constrain carbon uptake at the measured rate
        model_tmp = setParam(model_tmp,'lb',data.exchRxnIDs(exchIdx(i)),data.exchFluxes(i,exchIdx(i)));
    else % RMSE from max growth: leave carbon uptake fully open
        model_tmp = setParam(model_tmp,'lb',data.exchRxnIDs(exchIdx(i)),-1000);
    end

    % This row's own carbon source is never itself asserted zero-flux.
    rxn2blockIter = rxn2block;
    rxn2blockIter(strcmp(rxn2blockIter,data.exchRxnIDs(exchIdx(i)))) = [];

    o2Flux = data.exchFluxes(i,strcmp('oxygen',data.exchMets));
    if o2Flux == 0
        model_tmp = modelAdapter.makeModelAnaerobic(model_tmp);
    end
    if ~(isnan(data.Ptot(i)) || data.Ptot(i) == 0)
        model_tmp = modelAdapter.changeProteinBiomass(model_tmp,data.Ptot(i));
    end

    sol = solveLP(model_tmp);
    if checkSolution(sol)
        bioRxn  = getIndexes(model_tmp,data.biomass,'rxns');
        % Biomass carries an assumed 41 Cmmol/gDCW, matching computeExcarbon.
        bioMeas = data.grRate(i) * 41;
        bioSim  = sol.x(bioRxn) * 41;

        if constrain
            fluxToCheck   = ~isnan(data.exchFluxes(i,:));
            fluxRxnIdx    = getIndexes(model_tmp,data.exchRxnIDs(fluxToCheck),'rxns');
            cNormMeasFlux = model_tmp.excarbon(fluxRxnIdx) .* transpose(data.exchFluxes(i,fluxToCheck));
            cNormSimFlux  = model_tmp.excarbon(fluxRxnIdx) .* sol.x(fluxRxnIdx);

            blockRxnIdx = getIndexes(model_tmp,rxn2blockIter,'rxns');
            blockSim    = sol.x(blockRxnIdx) .* model_tmp.excarbon(blockRxnIdx);
            % Already-correct zero predictions don't dilute the RMSE.
            blockSim(blockSim == 0) = [];
            blockMeas = zeros(numel(blockSim),1);

            measured  = [cNormMeasFlux; bioMeas; blockMeas];
            simulated = [cNormSimFlux; bioSim; blockSim];
        else
            measured  = bioMeas;
            simulated = bioSim;
        end
        rmseList(i) = sqrt(mean((measured - simulated).^2));
    else
        rmseList(i) = NaN;
    end
end
rmseList(isnan(rmseList)) = 99; % Fixed penalty for an infeasible condition.
if isfield(data,'evotuneRMSEweight') && ~isempty(data.evotuneRMSEweight)
    rmseList = rmseList .* data.evotuneRMSEweight;
end
rmse = mean(rmseList);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = computeExcarbon(model)
% Per-reaction carbon-weight lookup for RMSE weighting, ported from
% fillCarbonNum.m. The biomass reaction gets 41 Cmmol/gDCW; every other
% exchange reaction gets the carbon count of its exchanged metabolite. A
% missing/unparseable formula falls back to 1 (handled by the caller, via
% ecModel.excarbon(ecModel.excarbon == 0) = 1, matching MATLAB's
% Ematrix(isnan(Ematrix)) = 1 inside getElementalComposition).
model.excarbon = zeros(length(model.rxns),1);
[EXrxn, EXrxnIdx] = getExchangeRxns(model);
CarbonNum = getcarbonnum(model,EXrxn);
model.excarbon(EXrxnIdx) = CarbonNum;
model.excarbon(strcmp(model.rxnNames,'growth')) = 41;
end

function CarbonNum = getcarbonnum(model,exrxn)
% exrxn must be exchange reactions containing exactly one metabolite.
[~,idx] = ismember(exrxn,model.rxns);
EXmets = model.S(:,idx);
EXmetsIdx = zeros(length(exrxn),1);
for k = 1:length(EXmets(1,:))
    EXmetsIdx(k) = find(EXmets(:,k));
end
EXfors = model.metFormulas(EXmetsIdx);
Ematrix = getElementalComposition(EXfors,{'C'});
Ematrix = Ematrix(:,1);
Ematrix(isnan(Ematrix)) = 1;
CarbonNum = Ematrix;
end

function [Ematrix, elements] = getElementalComposition(formulae, elements, chargeInFormula)
% Get the complete elemental composition matrix. It supports formulae with
% generic elements, parentheses and decimal places
%
% USAGE:
%    [Ematrix, elements] = getElementalComposition(formulae, elements, chargeInFormula)
%
% INPUT:
%    formulae:        cell array of strings of chemical formulae. Can contain any generic elements starting
%                     with a capital letter followed by lowercase letters or '_', followed by a non-negative number.
%                     Also support '()', '[]', '{}'. E.g. {'H2O'; '[H2O]2(CuSO4)Generic_element0.5'}
% OPTIONAL INPUTS:
%    elements:        elements from previous call to preserve the order (default {})
%    chargeInFormula: true to accept formulae containing the generic element 'Charge' representing the charges,
%                     followed by a real number, e.g., 'HCharge1', 'SO4Charge-2' (default false).
%
% OUTPUTS:
%    Ematrix:         elemental composition matrix (#formulae x #elements)
%    elements:         cell array of elements corresponding to the columns of Ematrix
%
% E.g., [Ematrix, elements] = getElementalComposition({'H2O'; '[H2O]2(CuSO4)Generic_element0.5'}) would return:
%  elements = {'H', 'O', 'Cu', 'S', 'Generic_element'}
%  Ematrix = [ 2,   1,    0,   0,   0;
%              4,   6,    1,   1,   0.5]
%
% Siu Hung Joshua Chan May 2017
% Vendored here (rather than assumed available on path) since it was
% previously only defined nested inside the now-removed fillCarbonNum.m.

if nargin < 3 || isempty(chargeInFormula)
    chargeInFormula = false;
else
    chargeInFormula = logical(chargeInFormula);
end
% for recalling the original formula at the top level if there are parentheses
% in the formula leading to iterative calling
persistent formTopLv
persistent formCurLv
% for storing error message during iterative calling
persistent errMsg
persistent errMsgInThisLoop
persistent topLvJ
persistent selfCall
if isempty(selfCall)
    selfCall = 0;
end
if ~selfCall
    [formTopLv, errMsg] = deal('');
end

if ~isstruct(formulae)
    if ~iscell(formulae)
        % make sure it is a cell array of strings
        formulae = {formulae};
    end
else
    % also accept COBRA model as input
    if isfield(formulae, 'metFormulas')
        formulae = formulae.metFormulas;
    else
        error('The 1st input ''formulae'' should be a cell array of strings of formulae or a COBRA model with *.metFormulas')
    end
end

if nargin < 2 || isempty(elements)
    elements = {};
elseif numel(unique(elements)) < numel(elements)
    error('Repeated elements in the input ''elements'' array.')
else
    elements = elements(:)';  % make sure it is a row vector
end

% the Ematrix
Ematrix = zeros(numel(formulae), numel(elements));
% replace all brackets and braces by parentheses
formulae = regexprep(formulae, '[\[\{]', '\(');
formulae = regexprep(formulae, '[\]\}]', '\)');
digit = floor(log10(numel(formulae))) + 1;
for j = 1:numel(formulae)
    if ~selfCall
        % reset top level information for each formula
        [formTopLv, topLvJ, errMsgInThisLoop] = deal(formulae{j}, j, false);
    end
    formulae{j} = strtrim(formulae{j});
    if ~isempty(formulae{j})
        % get all outer parentheses
        parenthesis = [];
        stP = [];
        stPpos = [];
        lv = 0;
        k = 1;
        while k <= length(formulae{j})
            if strcmp(formulae{j}(k),'(')
                if lv == 0
                    pStart = k;
                end
                lv = lv + 1;
            elseif strcmp(formulae{j}(k),')')
                if lv == 1
                    parenthesis = [parenthesis; [pStart k]];
                    % right parenthesis, get the following stoichiometry if any
                    stPre = regexp(formulae{j}(k+1:end), '^(\+|\-)?\d*\.?\d*', 'match');
                    if isempty(stPre)
                        stP = [stP; 1];
                        stPpos = [stPpos; k k];
                    else
                        stCheck = str2double(stPre{1});
                        s = '';
                        if isnan(stCheck)
                            s = 'Invalid';
                        elseif stCheck < 0
                            s = 'Negative';
                        end
                        if ~isempty(s)  % error if not convertible to number or negative
                            f = [formulae{j}(parenthesis(end,1):parenthesis(end,2)), stPre{1}];
                            addErrorMessage(sprintf('    %s stoichiometry in ''%s''\n', s, f))
                        end
                        % stoichiometry
                        stP = [stP; stCheck];
                        % position of the stoichiometry in the text
                        stPpos = [stPpos; k + 1, k + length(stPre{1})];
                        k = k + length(stPre{1});
                    end
                end
                lv = lv - 1;
            end
            k = k + 1;
        end
        if isempty(parenthesis)
        % No parenthesis, parse the formula
            re = regexp(formulae{j}, '([A-Z][a-z_]*)((?:\+|\-)?\d*\.?\d*)', 'tokens');
            s = strjoin(cellfun(@(x) strjoin(x,''), re, 'UniformOutput', false), '');
            errorFlag = 0;
            if ~strcmp(s, formulae{j})  % check if the entire formula is retrieved
                errorFlag = 1;
            else
                % check the stoichiometry
                eleCheck = cellfun(@(x) x{1}, re, 'UniformOutput', false);
                stCheck = cellfun(@(x) str2double(x{2}), re);
                stCheck(cellfun(@(x) isempty(x{2}), re)) = 1;  % empty value => stoichiometry = 1
                if any(isnan(stCheck))
                    f = strjoin(cellfun(@(x) strjoin(x, ''), re(isnan(stCheck)), 'UniformOutput', false), ''', ''');
                    [errorFlag, errMsgKey] = deal(2, 'Invalid');
                else
                    % if chargeInFormula is true, detect negative stoich for non-charge elements, else any negative stoich
                    negSt = (~chargeInFormula | ~strcmp(eleCheck, 'Charge')) & stCheck < 0;
                    if any(negSt)
                        f = strjoin(cellfun(@(x) strjoin(x, ''), re(negSt), 'UniformOutput', false), ''', ''');
                        [errorFlag, errMsgKey] = deal(2, 'Negative');
                    end
                end
            end
            if errorFlag > 0
                if selfCall <= 1
                    s2 = sprintf('from the part ''%s''', formulae{j});
                elseif selfCall == 2
                    s2 = sprintf('from the part ''(%s)''', formCurLv);
                elseif selfCall == 3
                    s2 = sprintf('from the part ''(%s)''', formulae{j});
                end
                if errorFlag == 1
                    s2 = sprintf(['    Only ''%s'' can be recognized %s.\n'...
                        '       Each element should start with a capital letter followed by lower case letters'...
                        ' or ''_'' with indefinite length and followed by a number.\n'], s, s2);
                elseif errorFlag == 2
                    s2 = sprintf('    %s stoichiometry in ''%s'' %s\n', errMsgKey, f, s2);
                end
                addErrorMessage(s2);
            end
            elementJ = repmat({''}, 1, numel(re));
            nEj = 0;
            stoichJ = zeros(numel(re),1);
            for k = 1:numel(re)
                [ynK,idK] = ismember(re{k}(1), elementJ(1:nEj));
                if ynK
                    k2 = idK;
                else
                    nEj = nEj + 1;
                    elementJ{nEj} = re{k}{1};
                    k2 = nEj;
                end
                if isempty(re{k}{2})
                    stoichJ(k2) = stoichJ(k2) + 1;
                else
                    stoichJ(k2) = stoichJ(k2) + str2double(re{k}{2});
                end
            end
            elementJ = elementJ(1:nEj);
            stoichJ = stoichJ(1:nEj);
            [ynE, idE] = ismember(elementJ, elements);  % map to existing elements
            if any(~ynE)
                idE(~ynE) = (numel(elements) + 1):(numel(elements) + sum(~ynE));
                Ematrix(:, (numel(elements) + 1):(numel(elements) + sum(~ynE))) = 0;
                elements = [elements, elementJ(~ynE)];
            end
            Ematrix(j, idE) = stoichJ;
        else
            % parentheses found. Iteratively get the formula inside parentheses
            rest = true(length(formulae{j}),1);
            for k = 1:size(parenthesis,1)
                selfCallCur = selfCall;
                selfCall = 3;
                [EmatrixK, elements] = getElementalComposition(formulae{j}(...
                    (parenthesis(k,1)+1):(parenthesis(k,2)-1)), elements, chargeInFormula);
                if numel(elements) > size(Ematrix, 2)
                    Ematrix(:, (size(Ematrix, 2) + 1):numel(elements)) = 0;
                end
                selfCall = selfCallCur;
                Ematrix(j, 1:numel(elements)) = Ematrix(j, 1:numel(elements)) + EmatrixK * stP(k);
                rest(parenthesis(k,1):stPpos(k,2)) = false;
            end
            if any(rest)
                formCurLv = formulae{j};
                selfCallCur = selfCall;
                selfCall = 1 + (selfCall >= 1);
                while any(rest)
                    % get the consecutive part of the formula that is not in parentheses
                    [pStart, pEnd] = deal(find(rest, 1), find(~rest));
                    pEnd = min(pEnd(pEnd > pStart));
                    if isempty(pEnd)
                        pEnd = numel(rest);
                    else
                        pEnd = pEnd - 1;
                    end
                    [EmatrixK, elements] = getElementalComposition(formulae{j}(pStart:pEnd), elements, chargeInFormula);
                    if numel(elements) > size(Ematrix, 2)
                        Ematrix(:, (size(Ematrix, 2) + 1):numel(elements)) = 0;
                    end
                    Ematrix(j, 1:numel(elements)) = Ematrix(j, 1:numel(elements)) + EmatrixK;
                    rest(pStart:pEnd) = false;
                end
                selfCall = selfCallCur;
            end
        end
    else
        Ematrix(j,:) = NaN;
    end
end
Ematrix = Ematrix(:, 1:numel(elements));
if ~selfCall
    selfCall = [];
    if ~isempty(errMsg)
        error(['%s\n', errMsg], 'Invalid formula input:')
    end
end

% nested function for adding error messages
    function addErrorMessage(s)
        if ~errMsgInThisLoop
            errMsg = [errMsg, sprintf(['#%0' num2str(digit) 'd:  %s\n'], topLvJ, formTopLv)];
            errMsgInThisLoop = true;
        end
        errMsg = [errMsg, s];
    end
end
