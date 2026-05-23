function plotBayesianDiagnostics(rmseTrace, kcatTrace, sigmaLogTrace, diagnostics, ecModel, modelAdapter)
% plotBayesianDiagnostics
%   Focused visualization of Bayesian kcat tuning results, emphasizing
%   which parameters changed, when they changed, and their impact on RMSE.
%
% Input:
%   rmseTrace      [1 x nGen+1] best RMSE per generation (incl. prior at index 1)
%   kcatTrace      [D x nGen+1] posterior kcat means per generation
%   sigmaLogTrace  [D x nGen+1] posterior sigmaLog per generation (optional)
%   diagnostics    struct from bayesianSensitivityTuning
%   ecModel        final ecModel with best kcats
%   modelAdapter   model adapter (optional, for source labels)

if nargin < 6
    modelAdapter = ModelAdapterManager.getDefault();
end

%% Handle multi-round merging (same as before)
roundBoundaries = [];
if iscell(rmseTrace)
    nRounds = length(rmseTrace);
    if ~iscell(kcatTrace) || length(kcatTrace) ~= nRounds
        error('When passing multiple rounds, rmseTrace and kcatTrace must be cell arrays of the same length.');
    end

    mergedRmse = rmseTrace{1}(:)';
    mergedKcat = kcatTrace{1};
    
    for r = 2:nRounds
        roundBoundaries(end+1) = length(mergedRmse) - 1;
        mergedRmse = [mergedRmse, rmseTrace{r}(2:end)];
        mergedKcat = [mergedKcat, kcatTrace{r}(:, 2:end)];
    end
    
    rmseTrace = mergedRmse;
    kcatTrace = mergedKcat;
end

%% Extract key metrics
nGenTotal = length(rmseTrace);
nGen = nGenTotal - 1;
genAxis = 0:nGen;

D = size(kcatTrace, 1);
kcat0 = kcatTrace(:, 1);
kcatFinal = kcatTrace(:, end);
logChanges = log10(kcatFinal ./ kcat0);

% Get source labels
if isfield(diagnostics, 'sourceLabels')
    sourceLabels = diagnostics.sourceLabels;
else
    sourceLabels = {'OpenKineticsPredictor', 'brenda', 'custom', 'unlabelled'};
end
nSources = length(sourceLabels);
colors = lines(nSources);

% Get source indices
kcatSourceIdx = zeros(size(ecModel.ec.kcat));
if isfield(ecModel, 'ec') && isfield(ecModel.ec, 'source')
    for i = 1:nSources-1
        idx = strcmpi(regexprep(ecModel.ec.source, '\s*\(.*$', ''), sourceLabels{i});
        kcatSourceIdx(idx) = i;
    end
end

% Get protein/reaction IDs and names
if isfield(ecModel, 'ec') && isfield(ecModel.ec, 'rxns')
    paramIDs = ecModel.ec.rxns;
else
    paramIDs = cellstr(num2str((1:D)'));  % Fallback to indices
end

if isfield(ecModel, 'ec') && isfield(ecModel.ec, 'rxnNames')
    paramNames = ecModel.ec.rxnNames;
elseif isfield(ecModel, 'rxnNames')
    % Try to map ec.rxns to rxnNames
    if isfield(ecModel, 'ec') && isfield(ecModel.ec, 'rxns')
        paramNames = cell(size(paramIDs));
        for i = 1:length(paramIDs)
            rxnIdx = find(strcmp(ecModel.rxns, paramIDs{i}), 1);
            if ~isempty(rxnIdx)
                paramNames{i} = ecModel.rxnNames{rxnIdx};
            else
                paramNames{i} = '';
            end
        end
    else
        paramNames = cell(size(paramIDs));
        paramNames(:) = {''};
    end
else
    paramNames = cell(size(paramIDs));
    paramNames(:) = {''};
end

% Define "changed" threshold (2% change in linear space, ~0.01 in log10 space)
changeThreshold = 0.01;
changed = abs(logChanges) > changeThreshold;

%% FIGURE 1: RMSE Convergence & Impact Analysis (Combined)
fig1 = figure('Position', [100, 100, 1400, 900], 'Name', 'RMSE & Impact Analysis');

% 1A: RMSE trace
subplot(2, 3, 1);
hold on;
plot(genAxis, rmseTrace, 'b-', 'LineWidth', 2.5);
xlabel('Generation', 'FontSize', 11);
ylabel('RMSE', 'FontSize', 11);
title('RMSE Evolution', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
addRoundMarkers(roundBoundaries);
ylim([0, max(rmseTrace) * 1.1]);
hold off;

% 1B: Per-generation improvement
subplot(2, 3, 2);
hold on;
rmseImprovement = [0, -diff(rmseTrace) ./ rmseTrace(1:end-1) * 100];
bar(genAxis, rmseImprovement, 'FaceColor', [0.3 0.6 0.9]);
xlabel('Generation', 'FontSize', 11);
ylabel('RMSE Improvement (%)', 'FontSize', 11);
title('Per-Generation RMSE Improvement', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
yline(0, 'r--', 'LineWidth', 1);
addRoundMarkers(roundBoundaries);
hold off;

% 1C: Top 10 most changed parameters with IDs (by magnitude, not impact)
subplot(2, 3, 3);
[sortedChanges, sortIdx] = sort(abs(logChanges), 'descend');
top10Idx = sortIdx(1:min(10, length(sortIdx)));
top10Changes = logChanges(top10Idx);
top10Sources = kcatSourceIdx(top10Idx);
top10Sources(top10Sources == 0) = nSources;
top10IDs = paramIDs(top10Idx);

% Escape underscores in IDs
top10IDsEscaped = cell(size(top10IDs));
for i = 1:length(top10IDs)
    top10IDsEscaped{i} = strrep(top10IDs{i}, '_', '\_');
end

barColors = colors(top10Sources, :);
barh(1:length(top10Idx), top10Changes);
set(gca, 'YDir', 'reverse');

% Add parameter IDs as y-tick labels
if length(top10IDsEscaped) <= 10
    set(gca, 'YTick', 1:length(top10IDsEscaped), 'YTickLabel', top10IDsEscaped);
end

xlabel('log_{10}(final/prior)', 'FontSize', 10);
ylabel('Parameter ID', 'FontSize', 10);
title('Top 10 Largest Kcat Changes (by magnitude)', 'FontSize', 11, 'FontWeight', 'bold');
grid on;
xline(0, 'k--', 'LineWidth', 1);

% Color bars by source
h = get(gca, 'Children');
for i = 1:length(h)
    if isa(h(i), 'matlab.graphics.chart.primitive.Bar')
        h(i).FaceColor = 'flat';
        h(i).CData = barColors;
    end
end

% 1D-1F: When did parameters change most (one subplot per source, including unlabelled)
maxSourcePlots = min(3, nSources);
for plotIdx = 1:maxSourcePlots
    subplot(2, 3, 3 + plotIdx);
    
    % Prioritize showing unlabelled if it exists, otherwise show first sources
    if nSources >= 4 && plotIdx == 3
        sourceIdx = nSources;  % Show unlabelled as third plot
    else
        sourceIdx = plotIdx;
    end
    
    if sourceIdx == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == sourceIdx;
    end
    
    sourceChanged = changed & idx;
    
    if sum(sourceChanged) == 0
        text(0.5, 0.5, 'No changes', 'HorizontalAlignment', 'center');
        title(sourceLabels{sourceIdx}, 'FontSize', 11, 'FontWeight', 'bold');
        continue;
    end
    
    % Find generation of maximum change for each parameter
    genOfMaxChange = zeros(D, 1);
    for i = find(sourceChanged)'
        trajectory = log10(kcatTrace(i, :) / kcat0(i));
        [~, genOfMaxChange(i)] = max(abs(diff(trajectory)));
    end
    
    histogram(genOfMaxChange(sourceChanged), 0:5:nGen, ...
        'FaceColor', colors(sourceIdx, :), 'EdgeColor', 'none');
    xlabel('Generation', 'FontSize', 10);
    ylabel('Number of Parameters', 'FontSize', 10);
    title(sprintf('%s: When Changed (n=%d)', sourceLabels{sourceIdx}, sum(sourceChanged)), ...
        'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    addRoundMarkers(roundBoundaries);
end

%% FIGURE 2: Kcat Changes by Source
fig2 = figure('Position', [120, 80, 1400, 800], 'Name', 'Kcat Changes by Source');

% First pass: determine global y-limits for violin plots
allChangedValues = [];
for sourceIdx = 1:nSources
    if sourceIdx == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == sourceIdx;
    end
    
    sourceChanges = logChanges(idx);
    sourceChanged = changed(idx);
    
    if sum(sourceChanged) > 0
        allChangedValues = [allChangedValues; sourceChanges(sourceChanged)];
    end
end

if ~isempty(allChangedValues)
    globalYLim = [min(allChangedValues) * 1.2, max(allChangedValues) * 1.2];
else
    globalYLim = [-1, 1];
end

% Top row: Violin plots with individual points
for sourceIdx = 1:nSources
    subplot(2, nSources, sourceIdx);
    
    % Get indices for this source
    if sourceIdx == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == sourceIdx;
    end
    
    if sum(idx) == 0
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', 'FontSize', 10);
        title(sourceLabels{sourceIdx}, 'FontSize', 11, 'FontWeight', 'bold');
        continue;
    end
    
    % Get changes for this source (ONLY changed parameters)
    sourceChanges = logChanges(idx);
    sourceChanged = changed(idx);
    
    nTotal = sum(idx);
    nChanged = sum(sourceChanged);
    nUnchanged = nTotal - nChanged;
    
    % Violin plot + scatter of ONLY changed parameters
    if nChanged > 0
        changedValues = sourceChanges(sourceChanged);
        
        hold on;
        % Draw violin
        try
            [f, xi] = ksdensity(changedValues, 'BoundaryCorrection', 'reflection');
            f = f / max(f) * 0.3;  % Normalize width
            
            fill([1 - f, fliplr(1 + f)], [xi, fliplr(xi)], colors(sourceIdx, :), ...
                'FaceAlpha', 0.4, 'EdgeColor', 'none');
        catch
            % If density fails, skip violin
        end
        
        % Add individual points with jitter
        jitter = 0.1 * (rand(size(changedValues)) - 0.5);
        scatter(1 + jitter, changedValues, 15, colors(sourceIdx, :), ...
            'filled', 'MarkerFaceAlpha', 0.6);
        
        % Add median line (red)
        medVal = median(changedValues);
        plot([0.7, 1.3], [medVal, medVal], 'r-', 'LineWidth', 2.5);
        
        % Add zero reference line
        yline(0, 'k--', 'LineWidth', 1);
        
        ylabel('log_{10}(final/prior)', 'FontSize', 10);
        title(sourceLabels{sourceIdx}, 'FontSize', 11, 'FontWeight', 'bold');
        set(gca, 'XTick', []);
        ylim(globalYLim);
        grid on;
        
        % Add text annotation with counts
        ylims = ylim;
        text(1.4, ylims(2)*0.95, sprintf('Changed: %d (%.1f%%)', nChanged, 100*nChanged/nTotal), ...
            'FontSize', 9, 'HorizontalAlignment', 'left');
        text(1.4, ylims(2)*0.85, sprintf('Unchanged: %d', nUnchanged), ...
            'FontSize', 9, 'HorizontalAlignment', 'left', 'Color', [0.5 0.5 0.5]);
        
        hold off;
    else
        text(0.5, 0.5, sprintf('No changes\n(n=%d)', nTotal), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        title(sourceLabels{sourceIdx}, 'FontSize', 11, 'FontWeight', 'bold');
        ylim(globalYLim);
    end
end

% Bottom row: Summary statistics per source
for sourceIdx = 1:nSources
    subplot(2, nSources, nSources + sourceIdx);
    
    if sourceIdx == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == sourceIdx;
    end
    
    if sum(idx) == 0
        continue;
    end
    
    sourceChanges = logChanges(idx);
    sourceChanged = changed(idx);
    
    nTotal = sum(idx);
    nChanged = sum(sourceChanged);
    nIncreased = sum(sourceChanges > changeThreshold);
    nDecreased = sum(sourceChanges < -changeThreshold);
    
    % Calculate fold changes (only for changed parameters)
    if nChanged > 0
        changedValues = sourceChanges(sourceChanged);
        
        if nIncreased > 0
            increasedValues = sourceChanges(sourceChanges > changeThreshold);
            meanIncFold = mean(10.^increasedValues);
            medIncFold = median(10.^increasedValues);
        else
            meanIncFold = NaN;
            medIncFold = NaN;
        end
        
        if nDecreased > 0
            decreasedValues = sourceChanges(sourceChanges < -changeThreshold);
            meanDecFold = mean(10.^abs(decreasedValues));
            medDecFold = median(10.^abs(decreasedValues));
        else
            meanDecFold = NaN;
            medDecFold = NaN;
        end
        
        % Display as text table
        axis off;
        text(0.05, 0.95, sprintf('\\bfTotal:\\rm %d', nTotal), ...
            'FontSize', 10, 'VerticalAlignment', 'top', 'Units', 'normalized');
        text(0.05, 0.80, sprintf('\\bfChanged:\\rm %d (%.1f%%)', nChanged, 100*nChanged/nTotal), ...
            'FontSize', 10, 'VerticalAlignment', 'top', 'Units', 'normalized');
        
        if nIncreased > 0
            text(0.05, 0.60, sprintf('\\bfIncreased:\\rm %d', nIncreased), ...
                'FontSize', 9, 'VerticalAlignment', 'top', 'Units', 'normalized', 'Color', [0.8 0 0]);
            text(0.05, 0.48, sprintf('  Mean: %.1fx', meanIncFold), ...
                'FontSize', 9, 'VerticalAlignment', 'top', 'Units', 'normalized', 'Color', [0.8 0 0]);
            text(0.05, 0.36, sprintf('  Median: %.1fx', medIncFold), ...
                'FontSize', 9, 'VerticalAlignment', 'top', 'Units', 'normalized', 'Color', [0.8 0 0]);
        end
        
        if nDecreased > 0
            text(0.05, 0.20, sprintf('\\bfDecreased:\\rm %d', nDecreased), ...
                'FontSize', 9, 'VerticalAlignment', 'top', 'Units', 'normalized', 'Color', [0 0 0.8]);
            text(0.05, 0.08, sprintf('  Mean: %.1fx', meanDecFold), ...
                'FontSize', 9, 'VerticalAlignment', 'top', 'Units', 'normalized', 'Color', [0 0 0.8]);
            text(0.05, -0.04, sprintf('  Median: %.1fx', medDecFold), ...
                'FontSize', 9, 'VerticalAlignment', 'top', 'Units', 'normalized', 'Color', [0 0 0.8]);
        end
    else
        axis off;
        text(0.5, 0.5, sprintf('No changes\n(n=%d)', nTotal), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
end

%% FIGURE 3: Temporal Evolution of Kcats
fig3 = figure('Position', [140, 60, 1400, 900], 'Name', 'Kcat Evolution Over Time');

% Select generations to plot
genSnapshots = unique([1, round(nGen*0.25), round(nGen*0.5), round(nGen*0.75), nGen]);
genSnapshots = genSnapshots(genSnapshots >= 1 & genSnapshots <= nGen);
if ~isempty(roundBoundaries)
    genSnapshots = unique([genSnapshots, roundBoundaries]);
end
genSnapshots = sort(genSnapshots);
nSnaps = length(genSnapshots);

% One subplot per source showing density evolution
for sourceIdx = 1:min(nSources, 6)
    subplot(2, 3, sourceIdx);
    
    if sourceIdx == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == sourceIdx;
    end
    
    if sum(idx) == 0
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
        title(sourceLabels{sourceIdx}, 'FontSize', 11, 'FontWeight', 'bold');
        continue;
    end
    
    hold on;
    genColors = parula(nSnaps);
    
    for snapIdx = 1:nSnaps
        gen = genSnapshots(snapIdx);
        kcatSnap = kcatTrace(idx, gen+1);
        
        try
            [f, xi] = ksdensity(log10(kcatSnap), 'BoundaryCorrection', 'reflection');
            
            if gen == nGen
                lineStyle = '-';
                lineWidth = 2.5;
                displayName = sprintf('Gen %d (final)', gen);
            elseif any(gen == roundBoundaries)
                lineStyle = '--';
                lineWidth = 2.5;
                displayName = sprintf('Gen %d (round boundary)', gen);
            else
                lineStyle = '-';
                lineWidth = 1.5;
                displayName = sprintf('Gen %d', gen);
            end
            
            plot(xi, f, lineStyle, 'LineWidth', lineWidth, ...
                'Color', genColors(snapIdx, :), 'DisplayName', displayName);
        catch
            % Skip if density estimation fails
        end
    end
    
    xlabel('log_{10}(kcat) [s^{-1}]', 'FontSize', 10);
    ylabel('Density', 'FontSize', 10);
    title(sprintf('%s (n=%d, %.1f%% changed)', ...
        sourceLabels{sourceIdx}, sum(idx), 100*sum(changed(idx))/sum(idx)), ...
        'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    xlim([-2, 5]);
    hold off;
end

%% FIGURE 4: Early Changes (First 5 Generations)
fig4 = figure('Position', [160, 40, 1400, 900], 'Name', 'Early Parameter Changes (Gen 1-5)');

% For each of first 5 generations, find parameters that changed
earlyGens = 1:min(5, nGen);
allEarlyChanges = [];

for genIdx = 1:length(earlyGens)
    gen = earlyGens(genIdx);
    
    % Find parameters that changed significantly in this generation
    if gen == 1
        prevKcat = kcat0;
    else
        prevKcat = kcatTrace(:, gen);
    end
    currentKcat = kcatTrace(:, gen+1);
    
    genChange = log10(currentKcat ./ prevKcat);
    genChanged = abs(genChange) > changeThreshold;
    
    if sum(genChanged) > 0
        changedIdx = find(genChanged);
        for i = 1:length(changedIdx)
            idx = changedIdx(i);
            
            % Determine source
            if kcatSourceIdx(idx) == 0
                srcIdx = nSources;
            else
                srcIdx = kcatSourceIdx(idx);
            end
            
            allEarlyChanges = [allEarlyChanges; struct(...
                'generation', gen, ...
                'paramIdx', idx, ...
                'paramID', paramIDs{idx}, ...
                'paramName', paramNames{idx}, ...
                'source', sourceLabels{srcIdx}, ...
                'sourceIdx', srcIdx, ...
                'change', genChange(idx), ...
                'foldChange', 10^abs(genChange(idx)), ...
                'direction', sign(genChange(idx)), ...
                'kcatInitial', prevKcat(idx), ...
                'kcatFinal', currentKcat(idx))];
        end
    end
end

if ~isempty(allEarlyChanges)
    % Group reactions that differ only by _EXP_\d+ and have similar kcats
    groupedChanges = groupSimilarReactions(allEarlyChanges, 0.05);  % 5% tolerance
    
    % Sort by magnitude of change
    [~, sortIdx] = sort(abs([groupedChanges.change]), 'descend');
    topN = min(20, length(sortIdx));
    topChanges = groupedChanges(sortIdx(1:topN));
    
    % Create table visualization
    subplot(1, 1, 1);
    axis off;
    
    % Header
    headerY = 0.98;
    text(0.02, headerY, 'Gen', 'FontSize', 9, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.08, headerY, 'Source', 'FontSize', 9, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.20, headerY, 'Reaction ID', 'FontSize', 9, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.42, headerY, 'Reaction Name', 'FontSize', 9, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.67, headerY, 'Initial kcat', 'FontSize', 9, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.77, headerY, 'Final kcat', 'FontSize', 9, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.87, headerY, 'Fold', 'FontSize', 9, 'FontWeight', 'bold', 'Units', 'normalized');
    
    % Rows
    rowHeight = 0.04;
    for i = 1:length(topChanges)
        yPos = headerY - (i * rowHeight);
        tc = topChanges(i);
        
        % Generation
        text(0.02, yPos, sprintf('%d', tc.generation), 'FontSize', 8, 'Units', 'normalized');
        
        % Source (colored)
        text(0.08, yPos, tc.source, 'FontSize', 8, 'Units', 'normalized', ...
            'Color', colors(tc.sourceIdx, :));
        
        % Reaction ID (escaped, possibly with count)
        paramStr = strrep(tc.paramID, '_', '\_');
        if tc.count > 1
            paramStr = sprintf('%s (n=%d)', paramStr, tc.count);
        end
        if length(paramStr) > 25
            paramStr = [paramStr(1:22), '...'];
        end
        text(0.20, yPos, paramStr, 'FontSize', 7, 'Units', 'normalized', 'Interpreter', 'tex');
        
        % Reaction Name (escaped)
        nameStr = strrep(tc.paramName, '_', '\_');
        if length(nameStr) > 30
            nameStr = [nameStr(1:27), '...'];
        end
        text(0.42, yPos, nameStr, 'FontSize', 7, 'Units', 'normalized', 'Interpreter', 'tex');
        
        % Initial kcat
        text(0.67, yPos, sprintf('%.2f', tc.kcatInitial), 'FontSize', 8, 'Units', 'normalized');
        
        % Final kcat
        text(0.77, yPos, sprintf('%.2f', tc.kcatFinal), 'FontSize', 8, 'Units', 'normalized');
        
        % Fold change (colored by direction)
        if tc.direction > 0
            changeColor = [0.8 0 0];  % Red for increase
        else
            changeColor = [0 0 0.8];  % Blue for decrease
        end
        text(0.87, yPos, sprintf('%.1fx', tc.foldChange), 'FontSize', 8, ...
            'Units', 'normalized', 'Color', changeColor);
    end
    
    title(sprintf('Top %d Parameter Changes in First 5 Generations (grouped by similar reactions)', length(topChanges)), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Units', 'normalized', 'Position', [0.5, 1.0, 0]);
else
    text(0.5, 0.5, 'No significant changes in first 5 generations', ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'Units', 'normalized');
end

%% Summary Report
fprintf('\n=== Bayesian Kcat Tuning Summary ===\n');
fprintf('Total parameters: %d\n', D);
fprintf('Total generations: %d\n', nGen);
fprintf('RMSE: %.4f → %.4f (%.1f%% reduction)\n', ...
    rmseTrace(1), rmseTrace(end), 100*(rmseTrace(1)-rmseTrace(end))/rmseTrace(1));
fprintf('\nPer-source summary:\n');
for i = 1:nSources
    if i == nSources
        idx = kcatSourceIdx == 0;
    else
        idx = kcatSourceIdx == i;
    end
    
    if sum(idx) == 0, continue; end
    
    nTotal = sum(idx);
    nChanged = sum(changed(idx));
    nInc = sum(logChanges(idx) > changeThreshold);
    nDec = sum(logChanges(idx) < -changeThreshold);
    
    fprintf('  %s: %d total, %d changed (%.1f%%), %d↑, %d↓\n', ...
        sourceLabels{i}, nTotal, nChanged, 100*nChanged/nTotal, nInc, nDec);
end

if ~isempty(allEarlyChanges)
    fprintf('\nEarly changes (first 5 generations): %d parameters\n', length(allEarlyChanges));
end

fprintf('\n');

end  % main function


%% --- Local helpers ----------------------------------------------------

function addRoundMarkers(boundaries)
if isempty(boundaries), return; end
for b = 1:length(boundaries)
    xline(boundaries(b), ':', sprintf('R%d', b+1), ...
        'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], ...
        'LabelOrientation', 'horizontal', ...
        'HandleVisibility', 'off');
end
end

function grouped = groupSimilarReactions(changes, tolerance)
% Group reactions that differ only by _EXP_\d+ suffix and have similar kcat values
%
% Input:
%   changes: struct array with fields paramID, kcatInitial, kcatFinal, etc.
%   tolerance: relative tolerance for considering kcats similar (e.g., 0.05 = 5%)
%
% Output:
%   grouped: struct array where similar reactions are combined with count field

if isempty(changes)
    grouped = changes;
    return;
end

% Extract base reaction names (remove _EXP_\d+ suffix)
nChanges = length(changes);
baseNames = cell(nChanges, 1);
for i = 1:nChanges
    rxnID = changes(i).paramID;
    % Remove _EXP_\d+ pattern (underscore, EXP, underscore, one or more digits)
    baseName = regexprep(rxnID, '_EXP_\d+$', '');
    baseNames{i} = baseName;
end

% Find unique base names
[uniqueBase, ~, groupIdx] = unique(baseNames);

% Group by base name and check if kcats are similar
grouped = [];
for i = 1:length(uniqueBase)
    groupMembers = find(groupIdx == i);
    
    if length(groupMembers) == 1
        % Single member, keep as is
        member = changes(groupMembers(1));
        member.count = 1;
        grouped = [grouped; member];
    else
        % Multiple members, check if kcats are similar
        kcatsInitial = [changes(groupMembers).kcatInitial];
        kcatsFinal = [changes(groupMembers).kcatFinal];
        
        % Check if all initial and final kcats are within tolerance
        meanInitial = mean(kcatsInitial);
        meanFinal = mean(kcatsFinal);
        
        initialSimilar = all(abs(kcatsInitial - meanInitial) <= tolerance * meanInitial);
        finalSimilar = all(abs(kcatsFinal - meanFinal) <= tolerance * meanFinal);
        
        if initialSimilar && finalSimilar
            % Combine into single entry
            representative = changes(groupMembers(1));
            representative.paramID = uniqueBase{i};  % Use base name
            representative.kcatInitial = meanInitial;
            representative.kcatFinal = meanFinal;
            representative.change = log10(meanFinal / meanInitial);
            representative.foldChange = 10^abs(representative.change);
            representative.count = length(groupMembers);
            grouped = [grouped; representative];
        else
            % Kcats differ significantly, keep separate
            for j = 1:length(groupMembers)
                member = changes(groupMembers(j));
                member.count = 1;
                grouped = [grouped; member];
            end
        end
    end
end
end
