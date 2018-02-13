%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch_CarbonSources_simulations
% 
% Ivan Domenzain. created:       2017-09-21
% Ivan Domenzain. Last modified: 2018-01-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flux_dist, conditions] = BatchSimulations(ecModel_batch,...
                                                                iteration)

figure
axis square
file_name = 'growthRates_data_carbonSources.txt';
fID       = fopen(file_name);
data      = textscan(fID,'%s %s %f','delimiter','\t');
efe       = fclose('all'); 
media     = [{'YEP'}, {'MAA'}, {'Min'}];
count      = 1;
conditions = [];
for i=1:length(media)
    media_indexes    = indexes_string(data{2},media{i},false);
    gRates_exp{1}{i} = data{1}(media_indexes);
    gRates_exp{2}{i} = data{3}(media_indexes);
    gRates_sim{i}    = [];
    SSres = 0;
    SStot = 0;
    for j=1:length(gRates_exp{2}{i})
        
        model = ecModel_batch;
        c_source        = strcat(gRates_exp{1}{i}{j},' exchange (reversible)');   
        [model,pos]     = changeMedia_batch(model,c_source,media{i});
        gR_pos          = find(strcmpi(model.rxnNames,'growth'));
        model.c         = zeros(size(model.c));
        model.c(gR_pos) = 1;
        
        solution           = solveLP(model);
        if solution.stat ~= 1
            disp(c_source)
            disp(solution.stat)
        end
        model.lb(gR_pos)   = 0.999*solution.x(gR_pos);
        model.ub(gR_pos)   = solution.x(gR_pos);
        solution           = solveLP(model,1);
        gRates_sim{i}  = [gRates_sim{i};solution.x(gR_pos)];
        res                 = abs((gRates_exp{2}{i}(j)-gRates_sim{i}(j))/...
            gRates_exp{2}{i}(j))*100;
        SSres               = SSres + res;
        
%         disp(gRates_exp{1}{i}{j})
%         disp(gRates_exp{2}{i}(j))
%         disp(gRates_sim{i}(j))
%         disp(res)
        switch gRates_exp{1}{i}{j}
                case 'D-fructose'
                    c_tag = 'Fru';
                case 'D-glucose'
                    c_tag = 'Glu';
                case 'sucrose'
                    c_tag = 'Suc';
                case 'maltose'
                   c_tag = 'Mal';
                case 'acetate'
                    c_tag = 'Ace';
                case 'D-galactose'
                    c_tag = 'Gal';
                case 'glycerol'
                    c_tag = 'Gly';
                case 'ethanol'
                    c_tag = 'Eth'; 
                case 'raffinose'
                    c_tag = 'Raf'; 
                case 'alpha,alpha-trehalose'
                    c_tag = 'Tre';
                case 'D-mannose'
                    c_tag = 'Man';

        end 
        flux_dist(:,count) = solution.x;
        count              = count+1;
        conditions         = [conditions;...
                              {char(strcat(media(i),string('_'),c_tag))}];
        switch i
              case 1
                  marker = 'd';
              case 2
                  marker = 's';
              case 3
                  marker = 'o';
        end
        
        text(gRates_exp{2}{i}(j),gRates_sim{i}(j),c_tag,'FontSize',14)
        hold on

    end
    plot(gRates_exp{2}{i},gRates_sim{i},marker,'MarkerSize',15)
    hold on
    title(['Max growth rate on different carbon sources #' num2str(iteration)],...
                                         'FontSize',30,'FontWeight','bold')
    ylabel('\mu_{max} predicted [h^{-1}]','FontSize',30,'FontWeight','bold');
    xlabel('\mu_{max} experimental [h^{-1}]','FontSize',30,'FontWeight','bold');
    RSq(i)       = (SSres)/length(gRates_exp{2}{i});
    legendStr(i) = strcat(media(i),' / e_{av}=',num2str(RSq(i)),'%');
end
x1 = linspace(0,1,1000);
plot(x1,x1)
hold on
legend(legendStr)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.
function matching = indexes_string(cell_array,str,flag)
    matching  = strfind(cell_array,str);
    if flag
       matching = find(~cellfun(@isempty,matching),1);
    else
        matching = find(~cellfun(@isempty,matching));
    end

end
