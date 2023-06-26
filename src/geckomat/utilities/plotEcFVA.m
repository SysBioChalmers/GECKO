function plotEcFVA(minFlux, maxFlux)
% plotEcFVA
%   Plots cumulative distribution functions of ecFVA results from one or
%   model ecModel(s) and/or GEM(s).
%
% Input:
%   minFlux     vector of minimum flux values, coming from ecFVA. If
%               multiple ecModels/GEMs are to be visualized, then each
%               column represents a separate model.
%   maxFlux     vector of maximum flux values, matching minFlux.

numMods = size(minFlux,2);

% Ignore zero flux reactions
for i=1:numMods
    zeroFlux = abs(minFlux(:,i)) < 1e-10 & abs(maxFlux(:,i)) < 1e-10;
    minFlux(zeroFlux,i) = NaN;
    maxFlux(zeroFlux,i) = NaN;
end

% Calculate flux ranges
fluxRange = maxFlux - minFlux;

% Plot all together
hold on
legendText = cell(1,numel(numMods));
for i=1:numMods
    cdfplot(fluxRange(:,i))
    legendText{i} = (['Model #' num2str(i) ' (median: ' num2str(median(fluxRange(:,i),'omitnan')) ')']);
end
set(gca, 'XScale', 'log', 'Xlim', [1e-7 1e4])
set(findall(gca, 'Type', 'Line'), 'LineWidth', 2)
legend(legendText,  'Location','northwest')
title('Flux variability (cumulative distribution)');
xlabel('Variability range (mmol/gDCWh)');
ylabel('Cumulative distribution');
