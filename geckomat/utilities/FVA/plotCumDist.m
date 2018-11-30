function [y_param, stats_param] = plotCumDist(E_parameter,legends,titlestr,median)
for i=1:length(E_parameter)
    
    if nargin<4
        median = false;
    end
      
    if median
        str{i} = horzcat(legends{i},' (',num2str(length(E_parameter{i})),' / ',num2str(median(E_parameter{i})), ')');
    else
        str{i} = horzcat(legends{i},' (',num2str(length(E_parameter{i})),')');
    end
    [y_param(i),stats_param(i)] = cdfplot(E_parameter{i});
    title(titlestr)
    ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
    xlabel('Variability range [mmol/gDw h]','FontSize',30,'FontWeight','bold');
    set(gca, 'XScale', 'log')
    hold on
end
legend(y_param,str);
end