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
fluxRanges = cell(3,1);
% Ignore zero flux reactions
for i=1:numMods
    zeroFlux = abs(minFlux(:,i)) < 1e-10 & abs(maxFlux(:,i)) < 1e-10;
    minFlux(zeroFlux,i) = NaN;
    maxFlux(zeroFlux,i) = NaN;
    fluxRange = maxFlux(:,i) - minFlux(:,i);
    fluxRange(isnan(fluxRange)) = [];
    fluxRanges{i} = fluxRange;
end

% Plot all together
hold on
legendText = cell(1,numel(numMods));
for i=1:numMods
    fluxRange = fluxRanges{i};
    cdfplot(fluxRange)
    legendText{i} = (['Model #' num2str(i) ' (median: ' num2str(median(fluxRange,'omitnan')) ')']);
end
set(gca, 'XScale', 'log', 'Xlim', [1e-7 1e4])
set(findall(gca, 'Type', 'Line'), 'LineWidth', 2)
legend(legendText,  'Location','northwest')
title('Flux variability (cumulative distribution)');
xlabel('Variability range (mmol/gDCWh)');
ylabel('Cumulative distribution');
hold off
end

function cdfplot(X)
% cdfplot(X) 
% displays a plot of the Empirical Cumulative Distribution Function 
% (CDF) of the input array X in the current figure. The empirical 
% CDF y=F(x) is defined as the proportion of X values less than or equal to x.
% If input X is a matrix, then cdfplot(X) parses it to the vector and 
% displays CDF of all values.
%
% EXAMPLE:
% figure;
% cdfplot(randn(1,100));
% hold on;
% cdfplot(-log(1-rand(1,100)));
% cdfplot(sqrt(randn(1,100).^2 + randn(1,100).^2))
% legend('Normal(0,1) CDF', 'Exponential(1) CDF', 'Rayleigh(1) CDF', 4)

% Version 1.0
% Alex Podgaetsky, September 2003
% alex@wavion.co.il
%
% Revisions:
%       Version 1.0 -   initial version

tmp = sort(reshape(X,prod(size(X)),1));
Xplot = reshape([tmp tmp].',2*length(tmp),1);

tmp = [1:length(X)].'/length(X);
Yplot = reshape([tmp tmp].',2*length(tmp),1);
Yplot = [0; Yplot(1:(end-1))];

figure(gcf);
hp = plot(Xplot, Yplot);

ColOrd = get(gca, 'ColorOrder'); 
ord = mod(length(get(gca,'Children')), size(ColOrd,1)); 
set(hp, 'Color', ColOrd((ord==0) + (ord>0)*ord, :));
if ~ishold
     xlabel('X', 'FontWeight','b','FontSize',12);
     ylabel('F(X)', 'FontWeight','b','FontSize',12);
     title('Empirical CDF', 'FontWeight','b','FontSize',12);
     grid on;
end
end