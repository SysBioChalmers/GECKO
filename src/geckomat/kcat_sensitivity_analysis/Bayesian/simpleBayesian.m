function [ecModel, rmseTrace, kcatTrace, sigmaLogTrace, diagnostics, posteriorSamples] = simpleBayesian(ecModel,kcatSigmaLog,modelAdapter)
% bayesianSensitivityTuning (SIMPLIFIED VERSION)
%   Performs ABC-SMC Bayesian parameter estimation to find kcat values that
%   match experimental flux data. This simplified version uses diagonal
%   proposals (no correlation learning) for easier understanding and
%   maintenance.
%
% SIMPLIFICATIONS vs FULL VERSION:
%   - No low-rank covariance learning (no SVD, no correlation detection)
%   - Simple independent proposals per parameter
%   - Fewer hyperparameters (removed: alpha, tau, cExpl, rMax, adaptFrac)
%   - May converge ~30-50% slower but easier to understand
%
% HOW IT WORKS:
%   1. INITIALIZATION: Starts with prior kcat values from databases
%   2. SAMPLING: Proposes new kcats independently per parameter in log-space
%   3. EVALUATION: Runs FBA and computes RMSE vs experimental data
%   4. SELECTION: Keeps best samples, discards poor ones
%   5. CONVERGENCE: Applies source-specific regularization and stops at plateau
%
% Input:
%   ecModel         ecModel with prior kcats
%   kcatSigmaLog    (optional) initial log-space std dev per parameter
%   modelAdapter    (optional) model adapter
%
% Output:
%   ecModel           ecModel with best kcats applied
%   rmseTrace         history of best RMSE per generation
%   kcatTrace         history of posterior kcat means per generation
%   sigmaLogTrace     history of posterior sigmaLog per generation
%   diagnostics       detailed diagnostics structure
%   posteriorSamples  final accepted samples (kcats, rmse)

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

%% Load hyperparameters
params = modelAdapter.params;

% INITIAL UNCERTAINTY (defines prior distributions)
sigma0logDefault     = params.bayesian.sigma0logDefault;
kcatSources          = params.bayesian.kcatSources;
sigma0logSource      = params.bayesian.sigma0logSource;

% REGULARIZATION (controls how parameters can change from prior)
shrinkThrDefault     = params.bayesian.shrinkThrDefault;
shrinkThrSource      = params.bayesian.shrinkThrSource;
varianceCapDefault   = params.bayesian.varianceCapDefault;
varianceCapSource    = params.bayesian.varianceCapSource;
forcePriorThrDefault = params.bayesian.forcePriorThrDefault;
forcePriorThrSource  = params.bayesian.forcePriorThrSource;
sparsityThr          = params.bayesian.sparsityThreshold;

% SAMPLING SCHEDULE
scheduleGenerations = params.bayesian.scheduleGenerations;
scheduleSamples     = params.bayesian.scheduleSamples;

% ACCEPTANCE CRITERIA
targetAccept        = params.bayesian.targetAccept;
minKeep             = params.bayesian.minKeep;
maxKeep             = params.bayesian.maxKeep;

% STOPPING CRITERIA
rmseThreshold       = params.bayesian.rmseThreshold;
maxGenerations      = params.bayesian.maxGenerations;
maxRMSEplateau      = params.bayesian.maxRMSEplateau;

% Get source indices
kcatSourceIdx = zeros(size(ecModel.ec.kcat));
for i = 1:numel(kcatSources)
    idx = strcmpi(regexprep(ecModel.ec.source, '\s*\(.*$', ''), kcatSources{i});
    kcatSourceIdx(idx) = i;
end
uniqKcatParams = find(kcatSourceIdx);
noKcatSource   = find(kcatSourceIdx == 0);

%% Construct initial kcat variances
if nargin < 2 || isempty(kcatSigmaLog)
    sigma0log = sigma0logDefault * ones(size(ecModel.ec.kcat));
    sigma0log(uniqKcatParams) = sigma0logSource(kcatSourceIdx(uniqKcatParams));
end

%% Load experimental data
bayData = loadBayesianData(modelAdapter);

% Add carbon counts for RMSE weighting
if ~isfield(ecModel,'excarbon')
    ecModel = addCarbonNum(ecModel);
    ecModel.excarbon(ecModel.excarbon == 0) = 1;
end

%% Initialize diagnostics
diagnostics = struct();
D = numel(ecModel.ec.kcat);

diagnostics.shrinkageTrace = zeros(D, maxGenerations);
diagnostics.acceptanceRateTrace = zeros(1, maxGenerations);
diagnostics.epsilonTrace = zeros(1, maxGenerations);
diagnostics.nSamplesTrace = zeros(1, maxGenerations);
diagnostics.nAcceptedTrace = zeros(1, maxGenerations);
diagnostics.activeBySource = zeros(numel(kcatSources)+1, maxGenerations);
diagnostics.nearPriorBySource = zeros(numel(kcatSources)+1, maxGenerations);
diagnostics.meanDeviationBySource = zeros(numel(kcatSources)+1, maxGenerations);
diagnostics.varianceRatioBySource = zeros(numel(kcatSources)+1, maxGenerations);
diagnostics.sparsityCountTrace = zeros(1, maxGenerations);
diagnostics.sparsityFractionTrace = zeros(1, maxGenerations);
diagnostics.diversityTrace = zeros(1, maxGenerations);
diagnostics.meanAcceptedRMSE = zeros(1, maxGenerations);
diagnostics.medianAcceptedRMSE = zeros(1, maxGenerations);
diagnostics.proposalWidthTrace = zeros(1, maxGenerations);

%% Initialize ABC-SMC
kcats = ecModel.ec.kcat;
kcat0 = kcats;

rmse = abc_max(ecModel, bayData, modelAdapter);
fprintf('RMSE with prior kcats: %.2f.\n', rmse)

rmseTop = rmse;
kcatTop = kcats;
rmseTrace = rmse;
kcatTrace = kcats;
sigmaLogTrace = sigma0log;

PB1 = ProgressBar2(maxGenerations, ['Run through ' num2str(maxGenerations) ' generations'], 'gui');

tmpModel = setParam(ecModel,'lb',modelAdapter.params.c_source,0);

%% Initialize proposal width (will be adapted)
sigmaProp_log = sigma0log;

%% ABC-SMC LOOP
generation = 1;
plateauCount = 0;

while rmse > rmseThreshold
    if generation <= maxGenerations
        fprintf('Running iteration %d of %d. ', generation, maxGenerations)

        schedulePos = find(generation >= scheduleGenerations, 1, 'last');
        N = scheduleSamples(schedulePos);
        newRmse = zeros(N,1);

        %% Draw new kcat samples
        if generation == 1
            % First generation: simple lognormal prior sampling
            muLog = log(kcats) - 0.5 .* (sigma0log .^ 2);
            muMat = repmat(muLog, 1, N);
            sigMat = repmat(sigma0log, 1, N);
            randomKcats = lognrnd(muMat, sigMat);
        else
            % Later generations: sample from accepted population
            nParents = size(kcatTop, 2);

            % Select parents with minimal duplication
            if N <= nParents
                parent_idx = randperm(nParents, N);
            else
                nReps = floor(N / nParents);
                nRemainder = N - nReps * nParents;
                base_idx = repmat(1:nParents, 1, nReps);
                extra_idx = randperm(nParents, nRemainder);
                parent_idx = [base_idx, extra_idx];
                parent_idx = parent_idx(randperm(N));
            end

            parents = kcatTop(:, parent_idx);
            
            % SIMPLIFIED PROPOSAL: Independent per-parameter in log-space
            randomKcats = proposeSimple(parents, sigmaProp_log, ecModel.ec.kcat);
        end

        % Guarantee positivity
        if any(randomKcats <= 0, 'all')
            warning('randomKcats contains non-positive values; enforcing floor.')
            randomKcats = max(randomKcats, realmin);
        end

        %% Evaluate RMSE for each proposal
        PB2 = ProgressBar2(N, ['Simulate ' num2str(N) ' models with random kcats'], 'gui');
        parfor j = 1:N
            ecModelIter = tmpModel;
            ecModelIter.ec.kcat = randomKcats(:,j);
            ecModelIter = applyKcatConstraints(ecModelIter);
            newRmse(j) = abc_max(ecModelIter, bayData, modelAdapter);
            count(PB2)
        end

        % Remove unsolvable models
        zeroRmse = newRmse == 0;
        newRmse(zeroRmse) = [];
        randomKcats(:,zeroRmse) = [];

        % Combine with previous accepted samples
        rmse = [newRmse; rmseTop];
        kcat = [randomKcats, kcatTop];

        %% Acceptance step: ABC-SMC epsilon thresholding
        epsilon = prctile(rmse, targetAccept);
        acc_idx = find(rmse <= epsilon);

        % Ensure at least minKeep and at most maxKeep survive
        minCount = max(1, floor(minKeep * numel(rmse)));
        maxCount = max(1, floor(maxKeep * numel(rmse)));
        if numel(acc_idx) < minCount
            [~, ord] = sort(rmse, 'ascend');
            acc_idx = ord(1:minCount);
            epsilon = rmse(acc_idx(end));
        elseif numel(acc_idx) > maxCount
            [~,ord] = sort(rmse(acc_idx), 'ascend');
            acc_idx = acc_idx(ord(1:maxCount));
        end

        rmseTop = rmse(acc_idx);
        kcatTop = kcat(:, acc_idx);

        %% Update posterior kcat and sigmaLog
        logKcatTop = log(kcatTop);
        muLog = mean(logKcatTop, 2);
        kcatSigmaLog = std(logKcatTop, 1, 2);

        % Compute deviation from prior
        devFromPrior = abs(muLog - log(kcat0)) ./ sigma0log;

        % Source-specific shrinkage
        shrinkWeight = min(devFromPrior / shrinkThrDefault, 1);
        for i = 1:numel(kcatSources)
            sourceIdx = kcatSourceIdx == i;
            shrinkWeight(sourceIdx) = min(devFromPrior(sourceIdx) / shrinkThrSource(i), 1);
        end

        % Force to prior if deviation too small
        forceToExactPrior = false(size(devFromPrior));
        for i = 1:numel(kcatSources)
            sourceIdx = (kcatSourceIdx == i);
            if any(sourceIdx) && forcePriorThrSource(i) > 0
                forceToExactPrior(sourceIdx) = (devFromPrior(sourceIdx) < forcePriorThrSource(i));
            end
        end
        if any(noKcatSource) && forcePriorThrDefault > 0
            forceToExactPrior(noKcatSource) = (devFromPrior(noKcatSource) < forcePriorThrDefault);
        end
        shrinkWeight(forceToExactPrior) = 0;

        % Bayesian update
        kcatSigmaLog_raw = shrinkWeight .* kcatSigmaLog + (1 - shrinkWeight) .* sigma0log;
        kcats_raw = exp(shrinkWeight .* muLog + (1 - shrinkWeight) .* log(kcat0));

        % Sparsity enforcement
        smalldiff = abs(log(kcats_raw) - log(kcat0)) < sparsityThr * sigma0log;
        kcats = kcats_raw;
        kcats(smalldiff) = kcat0(smalldiff);
        kcatSigmaLog = kcatSigmaLog_raw;
        kcatSigmaLog(smalldiff) = sigma0log(smalldiff);

        % Cap variance growth by source
        for i = 1:numel(kcatSources)
            sourceIdx = kcatSourceIdx == i;
            kcatSigmaLog(sourceIdx) = min(kcatSigmaLog(sourceIdx), sigma0log(sourceIdx) * varianceCapSource(i));
        end
        kcatSigmaLog(noKcatSource) = min(kcatSigmaLog(noKcatSource), sigma0log(noKcatSource) * varianceCapDefault);

        %% Update proposal width (SIMPLIFIED - no low-rank learning)
        sigmaProp_log = updateProposalWidth(log(kcatTop), sigma0log);

        % Track progress
        [bestRMSE, bestIdx] = min(rmseTop);
        rmseTrace = [rmseTrace, bestRMSE];
        kcatTrace = [kcatTrace, kcats];
        sigmaLogTrace = [sigmaLogTrace, kcatSigmaLog];

        % Early stopping check
        if generation >= 2
            improvement = abs(rmseTrace(end) - rmseTrace(end-1)) / rmseTrace(end-1);
            if improvement < 0.01
                plateauCount = plateauCount + 1;
            else
                plateauCount = 0;
            end

            if plateauCount >= maxRMSEplateau
                fprintf('├─────────────────────────────────────────────────────────────────────────────\n');
                fprintf('│ EARLY STOP: RMSE converged (plateau for %d generations)\n', plateauCount);
                fprintf('│ Final generation: %d │ Final RMSE: %.2f\n', generation, rmseTrace(end));
                fprintf('└─────────────────────────────────────────────────────────────────────────────\n');
                break;
            end
        end

        % Use best sample as center for next generation
        tmpModel.ec.kcat = kcatTop(:, bestIdx);
        tmpModel = applyKcatConstraints(tmpModel);

        %% Store diagnostics
        diagnostics.shrinkageTrace(:, generation) = shrinkWeight;
        diagnostics.acceptanceRateTrace(generation) = numel(rmseTop) / numel(rmse);
        diagnostics.epsilonTrace(generation) = epsilon;
        diagnostics.nSamplesTrace(generation) = N;
        diagnostics.nAcceptedTrace(generation) = numel(rmseTop);
        diagnostics.meanAcceptedRMSE(generation) = mean(rmseTop);
        diagnostics.medianAcceptedRMSE(generation) = median(rmseTop);
        diagnostics.sparsityCountTrace(generation) = sum(smalldiff);
        diagnostics.sparsityFractionTrace(generation) = sum(smalldiff) / D;

        diagnostics = storeSourceMetrics(diagnostics, generation, shrinkWeight, ...
            kcatSigmaLog, kcats, kcat0, sigma0log, kcatSourceIdx, ...
            kcatSources, noKcatSource);

        logKcatRanges = max(log(kcatTop), [], 2) - min(log(kcatTop), [], 2);
        diagnostics.diversityTrace(generation) = mean(logKcatRanges);
        diagnostics.proposalWidthTrace(generation) = mean(sigmaProp_log);

        printGenerationSummary(generation, rmseTrace, shrinkWeight, ...
            rmseTop, rmse, sum(smalldiff), D, ecModel, kcatSigmaLog, ...
            sigma0log, kcatSources)
        
        generation = generation + 1;
        count(PB1)
    else
        fprintf('Halted due to reaching maximum generation limit.\n')
        break
    end
end

%% Finalize
diagnostics.finalGeneration = generation;
diagnostics.converged = rmseTrace(end) < rmseThreshold;
diagnostics.sourceLabels = [kcatSources, {'unlabelled'}];

posteriorSamples = struct();
posteriorSamples.kcats = kcatTop;
posteriorSamples.rmse = rmseTop;

[~, bestIdx] = min(rmseTop);
ecModel.ec.kcat = kcatTop(:, bestIdx);
ecModel = applyKcatConstraints(ecModel);
fprintf('Final RMSE: %.2f.\n', rmseTop(bestIdx))
end

%% ========================================================================
%% SIMPLIFIED HELPER FUNCTIONS (replaces proposeLowRankMixture + buildLowRankLogProposal)
%% ========================================================================

function proposals = proposeSimple(parents, sigmaProp_log, priorKcat)
% proposeSimple - Simple independent proposals in log-space
%
% Replaces the complex proposeLowRankMixture function with diagonal
% proposals (no correlation learning).
%
% Input:
%   parents       [D x Nprop] current kcat values to perturb
%   sigmaProp_log [D x 1] proposal width per parameter (log-space)
%   priorKcat     [D x 1] prior kcat values (for biological bounds)
%
% Output:
%   proposals     [D x Nprop] new proposed kcat values

    [D, Nprop] = size(parents);
    
    % Independent normal steps in log-space
    logParents = log(max(parents, realmin));
    logSteps = randn(D, Nprop) .* sigmaProp_log;
    proposals = exp(logParents + logSteps);
    proposals = max(proposals, realmin);
    
    % Apply biological bounds (same as original)
    for i = 1:D
        prior_kcat = priorKcat(i);
        if prior_kcat > 1e4
            % Exceptional enzyme (e.g., catalase)
            proposals(i, :) = max(min(proposals(i, :), 1e8), prior_kcat / 100);
        else
            % Typical enzyme
            proposals(i, :) = max(min(proposals(i, :), 1e4), 1e-2);
        end
    end
end

function sigmaProp_log = updateProposalWidth(logKcatTop, sigma0log)
% updateProposalWidth - Adapt proposal width from accepted samples
%
% Replaces buildLowRankLogProposal but only returns marginal widths
% (no correlation structure).
%
% Input:
%   logKcatTop  [D x Nacc] log of accepted kcat samples
%   sigma0log   [D x 1] baseline prior width
%
% Output:
%   sigmaProp_log [D x 1] proposal width per parameter

    % Use observed spread (standard deviation)
    stds_obs = std(logKcatTop, 1, 2);
    
    % Blend with prior (50/50) and enforce floor
    sigmaProp_log = 0.5 * stds_obs + 0.5 * sigma0log;
    sigmaProp_log = max(sigmaProp_log, 0.15 * sigma0log);
end

%% ========================================================================
%% UNCHANGED HELPER FUNCTIONS (copied from original)
%% ========================================================================

function diagnostics = storeSourceMetrics(diagnostics, gen, shrinkWeight, ...
    kcatSigmaLog, kcats, kcat0, sigma0log, kcatSourceIdx, ...
    kcatSources, noKcatSource)
    
    for i = 1:numel(kcatSources)
        idx = (kcatSourceIdx == i);
        if any(idx)
            diagnostics.activeBySource(i, gen) = sum(shrinkWeight(idx) > 0.3);
            diagnostics.nearPriorBySource(i, gen) = sum(abs(kcatSigmaLog(idx) - sigma0log(idx)) < 0.1);
            diagnostics.meanDeviationBySource(i, gen) = mean(abs(log(kcats(idx)) - log(kcat0(idx))) ./ sigma0log(idx));
            diagnostics.varianceRatioBySource(i, gen) = mean(kcatSigmaLog(idx) ./ sigma0log(idx));
        end
    end
    
    otherIdx = numel(kcatSources) + 1;
    diagnostics.activeBySource(otherIdx, gen) = sum(shrinkWeight(noKcatSource) > 0.3);
    diagnostics.nearPriorBySource(otherIdx, gen) = sum(abs(kcatSigmaLog(noKcatSource) - sigma0log(noKcatSource)) < 0.1);
    diagnostics.meanDeviationBySource(otherIdx, gen) = mean(abs(log(kcats(noKcatSource)) - log(kcat0(noKcatSource))) ./ sigma0log(noKcatSource));
    diagnostics.varianceRatioBySource(otherIdx, gen) = mean(kcatSigmaLog(noKcatSource) ./ sigma0log(noKcatSource));
end

function printGenerationSummary(generation, rmseTrace, ...
    shrinkWeight, rmseTop, rmse, n_forced_prior, D, ...
    ecModel, kcatSigmaLog, sigma0log, kcatSources)

    acceptRate = numel(rmseTop) / numel(rmse);
    numel_rmseTop = numel(rmseTop);
    numel_rmse = numel(rmse);

    fprintf('│ RMSE: %5.2f →%5.2f (Δ%4.1f%%) │ Accept: %3d/%3d (%2.0f%%) │ Active: %4d (%2.0f%%) │ Sparse: %4d', ...
        rmseTrace(max(1, end-1)), rmseTrace(end), ...
        100 * (rmseTrace(max(1, end-1)) - rmseTrace(end)) / rmseTrace(max(1, end-1)), ...
        numel_rmseTop, numel_rmse, 100 * acceptRate, ...
        sum(shrinkWeight > 0.3), 100 * sum(shrinkWeight > 0.3) / D, ...
        n_forced_prior);

    warnings = {};
    if acceptRate > 0.5, warnings{end+1} = '⚠️ACC↑'; end
    if acceptRate < 0.08, warnings{end+1} = '⚠️ACC↓'; end
    if generation >= 5
        improvement = (rmseTrace(end-4) - rmseTrace(end)) / rmseTrace(end-4);
        if improvement < 0.02, warnings{end+1} = '⚠️PLAT'; end
    end
    if ~isempty(warnings)
        fprintf(' │ %s', strjoin(warnings, ' '));
    end
    fprintf('\n');

    if mod(generation, 10) == 0
        fprintf('├─────────────────────────────────────────────────────────────────────────────\n');
        fprintf('│ SOURCE BREAKDOWN:\n');
        
        for i = 1:numel(kcatSources)+1
            if i > numel(kcatSources)
                idx = ~ismember(regexprep(ecModel.ec.source, '\s*\(.*$', ''), kcatSources);
                sourceName = 'OTHERS';
            else
                idx = strcmpi(regexprep(ecModel.ec.source, '\s*\(.*$', ''), kcatSources{i});
                sourceName = upper(kcatSources{i});
            end

            if any(idx)
                n_total = sum(idx);
                n_active = sum(shrinkWeight(idx) > 0.3);
                n_near_prior = sum(abs(kcatSigmaLog(idx) - sigma0log(idx)) < 0.1);
                mean_var = mean(kcatSigmaLog(idx) ./ sigma0log(idx));
                
                fprintf('│   %-8s: %4d active (%2.0f%%) │ %4d near prior (%2.0f%%) │ σ: %.1fx  ', ...
                    sourceName, n_active, 100*n_active/n_total, ...
                    n_near_prior, 100*n_near_prior/n_total, mean_var);
                
                if (i == 1 && n_near_prior/n_total > 0.5) || ...
                   (i == 2 && n_active/n_total > 0.3) || ...
                   (i == 3 && n_near_prior/n_total > 0.7)
                    fprintf('✓\n');
                else
                    fprintf('⚠️\n');
                end
            end
        end
        
        if generation >= 10
            imp10 = 100 * (rmseTrace(end-9) - rmseTrace(end)) / rmseTrace(end-9);
            fprintf('│ 10-gen improvement: %.1f%%\n', imp10);
        end
        fprintf('└─────────────────────────────────────────────────────────────────────────────\n');
    end
end
