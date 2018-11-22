%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function model = block_production(model,metName,irrevFlag)
%  
% Function that blocks the secretion of the indicated metabolites in a GEM
% for avoiding weird solutions.
%
% INPUTS:
%    - model        A GEM matlab structure
%    - irrevFlag    True if the model is in a irreversible format
% OUTPUTS:
%    - model        The model structure with the modified constrains.
%
% Ivan Domenzain.      Last edited: 2018-03-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = block_production(model,metName,irrevFlag)
	% COBRA function for the extraction of Exchange reaction indexes
    [~,Uptake_indxs] = findExcRxns(model);
    Uptake_indxs 	 = find(Uptake_indxs);
    % Gets the indexes of all the extracellular mets
    compIndex		 = find(strcmpi(model.compNames,'extracellular'));
    compMets		 = find(model.metComps==compIndex);
    
    %For each metabolite uptake that is going to be blocked
    for i=1:length(metName)
        % Gets the index for the metabolite locallized in the extracellular
        % compartment
        ExtMet_Index = compMets(strcmpi(model.metNames(compMets),metName(i)));
        
        if length(ExtMet_Index)>1
            disp(['warning: multiple matches were found in the compartment for ' string(metName(i))])
        else
            
            % Get uptake rxns indexes for the especified metabolite
            metRxns      = find(model.S(ExtMet_Index,:));
            if irrevFlag==false
                metUptakeRxn = intersect(metRxns,Uptake_indxs);
                % Block metabolite secretion
                if ~isempty(metUptakeRxn)
                    model.ub(metUptakeRxn) = 0;
                end
            else
                Uptake_indxs = find(model.S(ExtMet_Index,:)<0);
                if ~isempty(Uptake_indxs)
                    for j=1:length(Uptake_indxs)
                        product = model.metNames(find(model.S(:,Uptake_indxs(j))>0));
                        if isempty(product)
                            break
                        end
                    end
                model.ub(Uptake_indxs(j)) = 0;
                end                
            end
        end
    end
end