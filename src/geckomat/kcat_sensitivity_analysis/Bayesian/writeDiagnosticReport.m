function writeDiagnosticReport(diagnostics, rmseTrace, kcatTrace, ecModel, modelAdapter, outputPath)
% writeDiagnosticReport - Generate comprehensive markdown report for Claude analysis
%
% This function creates a detailed markdown file containing all diagnostic
% information from a Bayesian kcat tuning run. The report can be fed to
% Claude in a new chat session to get advice on parameter tuning.
%
% Input:
%   diagnostics   - diagnostics struct from bayesianSensitivityTuning
%   rmseTrace     - RMSE trace from tuning
%   kcatTrace     - kcat trace from tuning
%   ecModel       - final ecModel with tuned kcats
%   modelAdapter  - model adapter with current parameter settings
%   outputPath    - (optional) path for output file
%                   Default: '/mnt/user-data/outputs/bayesian_diagnostic_report.md'
%
% Output:
%   Writes a markdown file ready to be pasted into Claude chat
%
% Example:
%   writeDiagnosticReport(diagnostics, rmseTrace, kcatTrace, ecModel, modelAdapter);
if nargin < 5 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
if nargin < 6 || isempty(outputPath)
    outputPath = fullfile(params.path,'output','bayesian_diagnostic_report.md');
end

%% Extract key information
nGen = diagnostics.finalGeneration - 1;
sourceLabels = diagnostics.sourceLabels;
nSources = length(sourceLabels);

%% Backward/forward compatibility
% Ensure every per-generation field this report consumes exists. Older runs
% may lack the newer traces (change-magnitude quartiles, new-proposal
% acceptance, boundary violations, directedness); fill any missing field
% with NaN so the report degrades gracefully instead of erroring.
nCol = max(nGen, 1);
twoDfields = {'changeMagnitudeQ25','changeMagnitudeQ50','changeMagnitudeQ75'};
for k = 1:numel(twoDfields)
    if ~isfield(diagnostics, twoDfields{k}) || isempty(diagnostics.(twoDfields{k}))
        diagnostics.(twoDfields{k}) = nan(nSources, nCol);
    end
end
oneDfields = {'proposalAcceptRate','boundaryViolationsTrace','directednessTrace'};
for k = 1:numel(oneDfields)
    if ~isfield(diagnostics, oneDfields{k}) || isempty(diagnostics.(oneDfields{k}))
        diagnostics.(oneDfields{k}) = nan(1, nCol);
    end
end

% Get source indices
kcatSourceIdx = zeros(size(ecModel.ec.kcat));
for i = 1:nSources-1
    idx = strcmpi(ecModel.ec.source, params.bayesian.kcatSources{i});
    kcatSourceIdx(idx) = i;
end

%% Open file for writing
fid = fopen(outputPath, 'w');
if fid == -1
    error('Could not open file: %s', outputPath);
end

%% Write header and instructions
fprintf(fid, '# Bayesian Kcat Tuning Diagnostic Report\n\n');
fprintf(fid, '**Generated:** %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '---\n\n');
fprintf(fid, '## Instructions for Claude\n\n');
fprintf(fid, 'Please analyze this diagnostic report from a Bayesian parameter tuning run and provide specific recommendations for improving the hyperparameter settings.\n\n');
fprintf(fid, 'Focus on:\n');
fprintf(fid, '1. Whether convergence was healthy or problematic\n');
fprintf(fid, '2. Whether source-specific regularization is working correctly\n');
fprintf(fid, '3. Specific hyperparameter adjustments with rationale\n');
fprintf(fid, '4. Any red flags or concerning patterns\n\n');
fprintf(fid, '---\n\n');

%% Current hyperparameter settings
fprintf(fid, '## Current Hyperparameter Settings\n\n');
fprintf(fid, '### Initial Uncertainty\n');
fprintf(fid, '```matlab\n');
fprintf(fid, 'sigma0logDefault = %.2f  %% Default log-space std dev\n', params.bayesian.sigma0logDefault);
fprintf(fid, 'kcatSources = {');
for i = 1:length(params.bayesian.kcatSources)
    fprintf(fid, '''%s''', params.bayesian.kcatSources{i});
    if i < length(params.bayesian.kcatSources), fprintf(fid, ', '); end
end
fprintf(fid, '}\n');
fprintf(fid, 'sigma0logSource = [');
for i = 1:length(params.bayesian.sigma0logSource)
    fprintf(fid, '%.2f', params.bayesian.sigma0logSource(i));
    if i < length(params.bayesian.sigma0logSource), fprintf(fid, '; '); end
end
fprintf(fid, ']  %% Per-source initial uncertainty\n');
fprintf(fid, '```\n\n');

fprintf(fid, '### Regularization (Source-Specific Constraints)\n');
fprintf(fid, '```matlab\n');
fprintf(fid, 'shrinkThrDefault = %.1f  %% Unlabelled: σ deviations for full update\n', params.bayesian.shrinkThrDefault);
fprintf(fid, 'shrinkThrSource = [');
for i = 1:length(params.bayesian.shrinkThrSource)
    fprintf(fid, '%.1f', params.bayesian.shrinkThrSource(i));
    if i < length(params.bayesian.shrinkThrSource), fprintf(fid, '; '); end
end
fprintf(fid, ']  %% Per-source shrinkage\n\n');

fprintf(fid, 'varianceCapDefault = %.2f  %% Unlabelled: max variance growth\n', params.bayesian.varianceCapDefault);
fprintf(fid, 'varianceCapSource = [');
for i = 1:length(params.bayesian.varianceCapSource)
    fprintf(fid, '%.2f', params.bayesian.varianceCapSource(i));
    if i < length(params.bayesian.varianceCapSource), fprintf(fid, '; '); end
end
fprintf(fid, ']  %% Per-source caps\n\n');

fprintf(fid, 'forcePriorThrDefault = %.1f  %% Unlabelled: snap to prior if dev < this\n', params.bayesian.forcePriorThrDefault);
fprintf(fid, 'forcePriorThrSource = [');
for i = 1:length(params.bayesian.forcePriorThrSource)
    fprintf(fid, '%.1f', params.bayesian.forcePriorThrSource(i));
    if i < length(params.bayesian.forcePriorThrSource), fprintf(fid, '; '); end
end
fprintf(fid, ']  %% Per-source force thresholds\n\n');

fprintf(fid, 'sparsityThreshold = %.2f  %% Clean up changes < this × sigma0\n', params.bayesian.sparsityThreshold);
fprintf(fid, '```\n\n');

fprintf(fid, '### Sampling & Acceptance\n');
fprintf(fid, '```matlab\n');
fprintf(fid, 'scheduleGenerations = [%s]\n', num2str(params.bayesian.scheduleGenerations));
fprintf(fid, 'scheduleSamples = [%s]\n', num2str(params.bayesian.scheduleSamples));
fprintf(fid, 'targetAccept = %.0f%%  %% Keep best X%% of samples\n', params.bayesian.targetAccept);
fprintf(fid, 'minKeep = %.2f  %% Min fraction to keep\n', params.bayesian.minKeep);
fprintf(fid, 'maxKeep = %.2f  %% Max fraction to keep\n', params.bayesian.maxKeep);
fprintf(fid, '```\n\n');

fprintf(fid, '### Stopping Criteria\n');
fprintf(fid, '```matlab\n');
fprintf(fid, 'rmseThreshold = %.2f  %% Target RMSE\n', params.bayesian.rmseThreshold);
fprintf(fid, 'maxGenerations = %d  %% Hard limit\n', params.bayesian.maxGenerations);
fprintf(fid, 'maxRMSEplateau = %d  %% Stop after this many generations without improvement\n', params.bayesian.maxRMSEplateau);
fprintf(fid, '```\n\n');

%% Run summary
fprintf(fid, '---\n\n');
fprintf(fid, '## Run Summary\n\n');
fprintf(fid, '- **Total generations:** %d\n', nGen);
fprintf(fid, '- **Converged:** %s\n', iif(diagnostics.converged, 'Yes', 'No'));
fprintf(fid, '- **Final RMSE:** %.4f\n', rmseTrace(end));
fprintf(fid, '- **Initial RMSE:** %.4f\n', rmseTrace(1));
fprintf(fid, '- **RMSE reduction:** %.1f%%\n', 100 * (rmseTrace(1) - rmseTrace(end)) / rmseTrace(1));
fprintf(fid, '- **Total parameters:** %d\n\n', numel(ecModel.ec.kcat));

%% Convergence analysis
fprintf(fid, '### Convergence Pattern\n\n');

% RMSE improvement by phase
earlyGens = 1:min(10, nGen);
midGens = max(11, min(10, nGen)+1):min(20, nGen);
lateGens = max(21, min(20, nGen)+1):nGen;

if length(earlyGens) > 1
    earlyImprovement = 100 * (rmseTrace(earlyGens(1)+1) - rmseTrace(earlyGens(end)+1)) / rmseTrace(earlyGens(1)+1);
    fprintf(fid, '- **Early phase (gen 1-10):** %.1f%% RMSE reduction\n', earlyImprovement);
end
if length(midGens) > 1
    midImprovement = 100 * (rmseTrace(midGens(1)+1) - rmseTrace(midGens(end)+1)) / rmseTrace(midGens(1)+1);
    fprintf(fid, '- **Mid phase (gen 11-20):** %.1f%% RMSE reduction\n', midImprovement);
end
if length(lateGens) > 1
    lateImprovement = 100 * (rmseTrace(lateGens(1)+1) - rmseTrace(lateGens(end)+1)) / rmseTrace(lateGens(1)+1);
    fprintf(fid, '- **Late phase (gen 21+):** %.1f%% RMSE reduction\n', lateImprovement);
end

% Direction consistency
if nGen >= 5
    lateDirectedness = mean(diagnostics.directednessTrace(max(1, end-9):end), 'omitnan');
    fprintf(fid, '- **Late-stage directedness:** %.2f (1=converging, 0.5=random, <0.3=oscillating)\n', lateDirectedness);
    if lateDirectedness < 0.3
        fprintf(fid, '  - ⚠️ **WARNING:** Strong oscillation detected\n');
    elseif lateDirectedness < 0.5
        fprintf(fid, '  - ⚠️ **WARNING:** Random walk behavior\n');
    end
end

fprintf(fid, '\n');

%% Per-source analysis
fprintf(fid, '---\n\n');
fprintf(fid, '## Per-Source Analysis\n\n');

for i = 1:nSources
    fprintf(fid, '### %s\n\n', sourceLabels{i});
    
    % Get indices for this source
    if i == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == i;
    end
    
    nTotal = sum(idx);
    if nTotal == 0
        fprintf(fid, '- No parameters from this source\n\n');
        continue;
    end
    
    % Final generation stats
    finalGen = nGen;
    nActive = diagnostics.activeBySource(i, finalGen);
    nNearPrior = diagnostics.nearPriorBySource(i, finalGen);
    meanDeviation = diagnostics.meanDeviationBySource(i, finalGen);
    varianceRatio = diagnostics.varianceRatioBySource(i, finalGen);
    
    % Change magnitudes
    q25 = diagnostics.changeMagnitudeQ25(i, finalGen);
    q50 = diagnostics.changeMagnitudeQ50(i, finalGen);
    q75 = diagnostics.changeMagnitudeQ75(i, finalGen);
    
    fprintf(fid, '**Summary:**\n');
    fprintf(fid, '- Total parameters: %d\n', nTotal);
    fprintf(fid, '- Active (shrinkWeight > 0.3): %d (%.1f%%)\n', nActive, 100*nActive/nTotal);
    fprintf(fid, '- Near prior (|σ_post - σ_prior| < 0.1): %d (%.1f%%)\n', nNearPrior, 100*nNearPrior/nTotal);
    fprintf(fid, '- Mean deviation from prior: %.2f σ\n', meanDeviation);
    fprintf(fid, '- Variance ratio (σ_post / σ_prior): %.2f\n', varianceRatio);
    fprintf(fid, '\n**Change magnitude (|log(final/prior)|):**\n');
    fprintf(fid, '- Q25: %.3f\n', q25);
    fprintf(fid, '- Median: %.3f\n', q50);
    fprintf(fid, '- Q75: %.3f\n\n', q75);
    
    % Expected behavior per source
    if i == nSources
        % Unlabelled - should be most active
        if nActive / nTotal < 0.20
            fprintf(fid, '⚠️ **WARNING:** Only %.0f%% active (expect 20-40%% for unlabelled)\n\n', 100*nActive/nTotal);
        end
    elseif strcmpi(sourceLabels{i}, 'brenda')
        % BRENDA - should stay mostly at prior
        if nActive / nTotal > 0.10
            fprintf(fid, '⚠️ **WARNING:** %.0f%% active (expect <5%% for BRENDA)\n\n', 100*nActive/nTotal);
        end
        if nNearPrior / nTotal < 0.90
            fprintf(fid, '⚠️ **WARNING:** Only %.0f%% near prior (expect >95%% for BRENDA)\n\n', 100*nNearPrior/nTotal);
        end
    elseif strcmpi(sourceLabels{i}, 'custom')
        % Custom - should be locked down
        if nNearPrior / nTotal < 0.95
            fprintf(fid, '⚠️ **WARNING:** Only %.0f%% near prior (expect ~100%% for custom)\n\n', 100*nNearPrior/nTotal);
        end
    end
    
    % Variance cap check
    if i <= length(params.bayesian.varianceCapSource)
        cap = params.bayesian.varianceCapSource(i);
        if abs(varianceRatio - cap) < 0.05
            fprintf(fid, '⚠️ **WARNING:** Variance ratio %.2f is at cap %.2f (may be constrained)\n\n', varianceRatio, cap);
        end
    elseif i == nSources
        cap = params.bayesian.varianceCapDefault;
        if abs(varianceRatio - cap) < 0.05
            fprintf(fid, '⚠️ **WARNING:** Variance ratio %.2f is at cap %.2f (may be constrained)\n\n', varianceRatio, cap);
        end
    end
end

%% Key diagnostic traces
fprintf(fid, '---\n\n');
fprintf(fid, '## Key Diagnostic Traces\n\n');

fprintf(fid, '### RMSE Evolution (every 5 generations)\n');
fprintf(fid, '```\n');
fprintf(fid, 'Gen | RMSE   | Improvement\n');
fprintf(fid, '----|--------|------------\n');
for g = [0, 5:5:nGen, nGen]
    if g == 0
        fprintf(fid, '%3d | %.4f | -\n', g, rmseTrace(g+1));
    else
        improvement = 100 * (rmseTrace(g) - rmseTrace(g+1)) / rmseTrace(g);
        fprintf(fid, '%3d | %.4f | %.1f%%\n', g, rmseTrace(g+1), improvement);
    end
end
fprintf(fid, '```\n\n');

fprintf(fid, '### Acceptance Rates (every 5 generations)\n');
fprintf(fid, '```\n');
fprintf(fid, 'Gen | Overall | New Proposals | Boundary Violations\n');
fprintf(fid, '----|---------|---------------|--------------------\n');
for g = 5:5:nGen
    fprintf(fid, '%3d | %5.1f%% | %5.1f%%       | %5.1f%%\n', ...
        g, ...
        100*diagnostics.acceptanceRateTrace(g), ...
        100*diagnostics.proposalAcceptRate(g), ...
        100*diagnostics.boundaryViolationsTrace(g));
end
fprintf(fid, '```\n\n');

fprintf(fid, '### Directedness & Diversity (every 5 generations)\n');
fprintf(fid, '```\n');
fprintf(fid, 'Gen | Directedness | Diversity | Sparsity %%\n');
fprintf(fid, '----|--------------|-----------|----------\n');
for g = 5:5:nGen
    fprintf(fid, '%3d | %4.2f        | %5.3f     | %5.1f%%\n', ...
        g, ...
        diagnostics.directednessTrace(g), ...
        diagnostics.diversityTrace(g), ...
        100*diagnostics.sparsityFractionTrace(g));
end
fprintf(fid, '```\n\n');

%% Problem indicators
fprintf(fid, '---\n\n');
fprintf(fid, '## Potential Issues Detected\n\n');

problemsFound = false;

% Check 1: Oscillation
if nGen >= 10
    lateDirectedness = mean(diagnostics.directednessTrace(max(1, end-9):end), 'omitnan');
    if lateDirectedness < 0.3
        fprintf(fid, '### ⚠️ CRITICAL: Strong Oscillation\n');
        fprintf(fid, '- Late-stage directedness: %.2f (should be >0.5)\n', lateDirectedness);
        fprintf(fid, '- **Likely cause:** Shrinkage thresholds too aggressive, parameters fighting between prior and posterior\n');
        fprintf(fid, '- **Recommendation:** Increase shrinkThrSource values by 20-30%%\n\n');
        problemsFound = true;
    elseif lateDirectedness < 0.5
        fprintf(fid, '### ⚠️ WARNING: Random Walk Detected\n');
        fprintf(fid, '- Late-stage directedness: %.2f (should be >0.5)\n', lateDirectedness);
        fprintf(fid, '- **Likely cause:** RMSE plateau reached, algorithm not stopping early enough\n');
        fprintf(fid, '- **Recommendation:** Decrease maxRMSEplateau from %d to %d\n\n', ...
            params.bayesian.maxRMSEplateau, max(5, params.bayesian.maxRMSEplateau - 3));
        problemsFound = true;
    end
end

% Check 2: High boundary violations
if nGen >= 5
    lateBoundary = mean(diagnostics.boundaryViolationsTrace(max(1, end-9):end));
    if lateBoundary > 0.1
        fprintf(fid, '### ⚠️ WARNING: High Boundary Violations\n');
        fprintf(fid, '- Late-stage violations: %.1f%% (should be <5%%)\n', 100*lateBoundary);
        fprintf(fid, '- **Likely cause:** Proposal width too large or variance caps too loose\n');
        fprintf(fid, '- **Recommendation:** Reduce varianceCapSource by 10-20%%\n\n');
        problemsFound = true;
    end
end

% Check 3: Poor acceptance rate
if nGen >= 5
    lateAccept = mean(diagnostics.proposalAcceptRate(max(1, end-9):end));
    if lateAccept < 0.03
        fprintf(fid, '### ⚠️ WARNING: Very Low Proposal Acceptance\n');
        fprintf(fid, '- Late-stage new proposal acceptance: %.1f%% (should be 5-15%%)\n', 100*lateAccept);
        fprintf(fid, '- **Likely cause:** Proposals too aggressive\n');
        fprintf(fid, '- **Recommendation:** Proposal width adaptation may need adjustment\n\n');
        problemsFound = true;
    elseif lateAccept > 0.3
        fprintf(fid, '### ⚠️ WARNING: Very High Proposal Acceptance\n');
        fprintf(fid, '- Late-stage new proposal acceptance: %.1f%% (should be 5-15%%)\n', 100*lateAccept);
        fprintf(fid, '- **Likely cause:** Proposals too conservative, not exploring enough\n');
        fprintf(fid, '- **Recommendation:** May need larger proposal width or fewer samples per generation\n\n');
        problemsFound = true;
    end
end

% Check 4: Source-specific issues
for i = 1:nSources
    if i == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == i;
    end
    
    nTotal = sum(idx);
    if nTotal == 0, continue; end
    
    finalGen = nGen;
    nActive = diagnostics.activeBySource(i, finalGen);
    nNearPrior = diagnostics.nearPriorBySource(i, finalGen);
    
    % BRENDA should stay locked
    if strcmpi(sourceLabels{i}, 'brenda') && nActive / nTotal > 0.10
        fprintf(fid, '### ⚠️ WARNING: BRENDA Too Active\n');
        fprintf(fid, '- Active BRENDA parameters: %.1f%% (should be <5%%)\n', 100*nActive/nTotal);
        fprintf(fid, '- **Recommendation:** Increase BRENDA shrinkThr from %.1f to %.1f\n\n', ...
            params.bayesian.shrinkThrSource(2), params.bayesian.shrinkThrSource(2) * 1.3);
        problemsFound = true;
    end
    
    % Custom should be completely locked
    if strcmpi(sourceLabels{i}, 'custom') && nNearPrior / nTotal < 0.98
        fprintf(fid, '### ⚠️ WARNING: Custom Values Not Locked\n');
        fprintf(fid, '- Custom parameters near prior: %.1f%% (should be ~100%%)\n', 100*nNearPrior/nTotal);
        fprintf(fid, '- **Recommendation:** Increase custom forcePriorThr from %.1f to %.1f\n\n', ...
            params.bayesian.forcePriorThrSource(3), params.bayesian.forcePriorThrSource(3) * 1.2);
        problemsFound = true;
    end
    
    % Unlabelled should be active
    if i == nSources && nActive / nTotal < 0.15
        fprintf(fid, '### ⚠️ WARNING: Unlabelled Parameters Too Constrained\n');
        fprintf(fid, '- Active unlabelled parameters: %.1f%% (should be 20-40%%)\n', 100*nActive/nTotal);
        fprintf(fid, '- **Recommendation:** Decrease shrinkThrDefault from %.1f to %.1f\n\n', ...
            params.bayesian.shrinkThrDefault, params.bayesian.shrinkThrDefault * 0.8);
        problemsFound = true;
    end
end

if ~problemsFound
    fprintf(fid, '✅ No major issues detected. Parameter settings appear reasonable.\n\n');
end

%% Detailed data section
fprintf(fid, '---\n\n');
fprintf(fid, '## Detailed Diagnostic Data\n\n');
fprintf(fid, '<details>\n');
fprintf(fid, '<summary>Click to expand full diagnostic traces</summary>\n\n');

% Full RMSE trace
fprintf(fid, '### Complete RMSE Trace\n');
fprintf(fid, '```\n');
for g = 0:nGen
    fprintf(fid, '%d,%.6f\n', g, rmseTrace(g+1));
end
fprintf(fid, '```\n\n');

% Acceptance rate trace
fprintf(fid, '### Acceptance Rate Trace\n');
fprintf(fid, '```\n');
fprintf(fid, 'Gen,Overall,NewProposals,BoundaryViolations\n');
for g = 1:nGen
    fprintf(fid, '%d,%.4f,%.4f,%.4f\n', ...
        g, ...
        diagnostics.acceptanceRateTrace(g), ...
        diagnostics.proposalAcceptRate(g), ...
        diagnostics.boundaryViolationsTrace(g));
end
fprintf(fid, '```\n\n');

% Directedness trace
fprintf(fid, '### Directedness Trace\n');
fprintf(fid, '```\n');
for g = 2:nGen
    fprintf(fid, '%d,%.4f\n', g, diagnostics.directednessTrace(g));
end
fprintf(fid, '```\n\n');

% Per-source change magnitudes
for i = 1:nSources
    fprintf(fid, '### %s Change Magnitude Quartiles\n', sourceLabels{i});
    fprintf(fid, '```\n');
    fprintf(fid, 'Gen,Q25,Q50,Q75\n');
    for g = 1:nGen
        fprintf(fid, '%d,%.4f,%.4f,%.4f\n', ...
            g, ...
            diagnostics.changeMagnitudeQ25(i, g), ...
            diagnostics.changeMagnitudeQ50(i, g), ...
            diagnostics.changeMagnitudeQ75(i, g));
    end
    fprintf(fid, '```\n\n');
end

fprintf(fid, '</details>\n\n');

%% Final prompt for Claude
fprintf(fid, '---\n\n');
fprintf(fid, '## Analysis Request\n\n');
fprintf(fid, 'Based on this diagnostic report, please provide:\n\n');
fprintf(fid, '1. **Overall assessment** of convergence quality\n');
fprintf(fid, '2. **Specific parameter recommendations** with justification:\n');
fprintf(fid, '   - Which shrinkThrSource values to adjust?\n');
fprintf(fid, '   - Which varianceCapSource values to adjust?\n');
fprintf(fid, '   - Should sparsityThreshold change?\n');
fprintf(fid, '   - Should maxRMSEplateau change?\n');
fprintf(fid, '3. **Priority order** for changes (what to fix first)\n');
fprintf(fid, '4. **Expected impact** of recommended changes\n\n');
fprintf(fid, 'Please be specific with numbers (e.g., "increase shrinkThrSource(2) from 10.0 to 12.0") rather than general advice.\n');

%% Close file
fclose(fid);

fprintf('Diagnostic report written to: %s\n', outputPath);
fprintf('Copy this file content into a Claude chat for tuning advice.\n');

end

function result = iif(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end
