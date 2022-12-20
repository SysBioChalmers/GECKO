%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [limKcats,limRxns,breakFlag] = TopLimitingKcat(model,prevUnicodes,...
%                                                                  prevRxns)
%
% Function that gets an EC model and returns the top limiting Kcat value in
% the network based on a sensitivity analysis on the objective function for 
% each of the assigned parameters (objective control coefficients)
%
% If the model is not feasible then the analysis is performed on a reaction 
% basis (relaxing all of its matched kcats) and gets the top limiting one, 
% if any.
%
% Ivan Domenzain    Last edited. 2018-09-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [limitations,breakFlag] = findTopLimitations(model,prevUnicodes,option)
                                                       
    %Get an initial simulation optimizing for biomass production
    base_sol = solveLP(model); 
    %Extract the metabolite indexes corresponding to the model enzymes
    %(avoiding the protein pool)
    enz_indxs = find(~cellfun(@isempty,strfind(model.metNames,'prot_')));
    enz_indxs = enz_indxs(1:end-1);
    %Extract the enzyme usage pseudoreactions indexes
    enzUsageIndxs = find(~cellfun(@isempty,strfind(model.rxns,'prot_'))); 
    enzUsageIndxs = enzUsageIndxs(1:end-1);
    limitations   = [];
    breakFlag     = false;
    %If the model is feasible, analyse the individual kcat values of the flux
    %carrier enzymes
    if ~isempty(base_sol.f)
        %First, find the limiting Kcats in the model
        limKcats = findLimitingKcat(model,base_sol,enzUsageIndxs);
        if ~isempty(limKcats{1})
            %Sort the limiting Kcats cell according to their control coefficient 
            limKcats = sortByCCoeff(limKcats,limKcats{5});
            j        = 1;
            pos      = 0;
            %Choose the top controlling Kcat that is not included in the cell
            %prevUnicodes to avoid infinite loops
            while j<=length(limKcats{1})
                %The top control coefficient is identified with its corresponding
                %uniprot code and the rxn index in the model because there might 
                %be promiscous enzymes catalysing different reactions at different rates
                Protname = [limKcats{1}{j}(6:end) '_' num2str(limKcats{3}(j))];
                if ~ismember(Protname,prevUnicodes)
                    pos = j;
                    break
                else
                    j = j+1;
                end
            end

            for i=1:length(limKcats)
                if pos~=0
                    limKcats{i} = limKcats{1,i}(pos);
                else
                    %If all the significant Limiting Kcats have been previously
                    %modified then the iterations should be stopped
                    limKcats{i} = [];
                    breakFlag   = true;
                end
            end
        end
        
        if ~isempty(limKcats{1})
            limitations = limKcats;
        else
            breakFlag = true;
        end
    %If the model is not feasible, analyse all of the kinetic parameters
    %associated to each of the metabolic reactions.
    else
        breakFlag = true;
        if option == 1
            limRxns     = findLimitingRxn(model,enz_indxs,enzUsageIndxs);
            limitations = limRxns;
        elseif option == 2
            limEnz      = findLimitingEnz(model,enz_indxs,enzUsageIndxs);
            limitations = limEnz;
        end
    end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_cell = sortByCCoeff(input_cell,values)
    [~,indexes] = sort(values,'descend');   
    for i=1:length(input_cell)
        output_cell{i} = input_cell{i}(indexes);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limKcats = findLimitingKcat(model,base_sol,enzUsageIndxs)
limKcats = cell(1,6);
%Get the original rxn indexes of the flux carrier enzymes (enzyme usages)
enzUsageFlux = base_sol.x(enzUsageIndxs);
fluxCarriers = find(enzUsageFlux);
fluxCarriers = enzUsageIndxs(fluxCarriers);
%Extract the metabolite index for the flux carrying enzymes
usedEnzIndxs = zeros(size(fluxCarriers));
for i=1:length(usedEnzIndxs)
    usedEnzIndxs(i) = find(model.S(:,fluxCarriers(i))==1);
end
%For each flux carrying enzyme calculate the objective control
%coefficient (ObjCC) on each of its kinetic coefficients
for i = 1:length(usedEnzIndxs)
    enzPos  = usedEnzIndxs(i);
    enzName = model.mets(enzPos);
    %Get the nonzero kinetic coefficients
    enz_rxnsCoeffs = find(model.S(enzPos,:));
    %For each enzymatic reaction (avoid enzyme usage pseudoRxn)
    if ~isempty(enz_rxnsCoeffs)
        for j=1:length(enz_rxnsCoeffs)-1
            temp_model = model;
            coeffPos   = enz_rxnsCoeffs(j);
            temp_model.S(enzPos,coeffPos) = model.S(enzPos,coeffPos)/1000;
            %Get a new solution with flexibilized j-th Kcat for the i-th
            %enzyme
            new_sol = solveLP(temp_model);
            if ~isempty(new_sol.f)
                %If solution was feasible then calculate the control
                %coefficient for the j-th Kcat
                objCC = (new_sol.f - base_sol.f)/base_sol.f;
                if abs(objCC) > 1e-3
                    Kcat         = (-1/model.S(enzPos,coeffPos))/3600;
                    limKcats{1} = [limKcats{1}; enzName];
                    limKcats{2} = [limKcats{2}; enzPos];
                    limKcats{3} = [limKcats{3}; coeffPos];
                    limKcats{4} = [limKcats{4}; Kcat];
                    limKcats{5} = [limKcats{5}; objCC];
                    limKcats{6} = [limKcats{6}; model.rxnNames(coeffPos)];
                end
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [limRxns] = findLimitingRxn(model,enz_indxs,enzUsageIndxs) 
limRxns   = cell(1,3);
for i=1:(enzUsageIndxs(1)-1)
    %disp(['Relaxing Kcat coefficients for Rxn#' num2str(i)])
    temp_model = model;
    %Search if there are nonZero Kcat coefficients matched to the
    %i-th reaction
    nonZero = find(temp_model.S(enz_indxs,i));
    if ~isempty(nonZero)
        disp(['Analyzing rxn: #' num2str(i)])
        %Flexibilize all non-zero coefficients in the i-th metabolic reaction 
        %and check if a non-zero value for the objective function can be obtained
        pos                 = enz_indxs(nonZero);
        temp_model.S(pos,i) = temp_model.S(pos,i)/1e6;
        new_sol             = solveLP(temp_model);
        new_sol             = abs(new_sol.f);
        if new_sol>0
            disp(['limiting rxn found: #' num2str(i)])
            disp(new_sol)
            limRxns{1} = [limRxns{1}; model.rxnNames(i)];
            limRxns{2} = [limRxns{2}; i];
            limRxns{3} = [limRxns{3}; new_sol];
        end
    end
end
%Sort results according to the predicted value for the objective function
if ~isempty(limRxns{3})
    limRxns = sortByCCoeff(limRxns,limRxns{3});
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [limEnz] = findLimitingEnz(model,enz_indxs,enzUsageIndxs)
limEnz   = cell(1,3);
for i=1:length(enz_indxs)
    enzPos = enz_indxs(i);
    %disp(['Relaxing Kcat coefficients for Enzyme #' num2str(i)])
    temp_model = model;
    %Search if there are nonZero Kcat coefficients matched to the
    %i-th enzyme
    nonZero = find(temp_model.S(enzPos,1:enzUsageIndxs(1)-1));
    if ~isempty(nonZero)
        disp(['Analyzing enzyme: #' num2str(i)])
        %Flexibilize all non-zero coefficients in the i-th enzyme and check 
        %if a non-zero value for the objective function can be obtained
        temp_model.S(enzPos,nonZero) = temp_model.S(enzPos,nonZero)/1e6;
        
        new_sol  = solveLP(temp_model);
        new_sol  = abs(new_sol.f);
        
        if new_sol>0
            disp(['limiting enzyme found: #' num2str(i)])
            disp(new_sol)
            limEnz{1} = [limEnz{1}; model.metNames(i)];
            limEnz{2} = [limEnz{2}; i];
            limEnz{3} = [limEnz{3}; new_sol];
        end
    end
end
%Sort results according to the predicted value for the objective function
if ~isempty(limEnz{3})
    limEnz = sortByCCoeff(limEnz,limEnz{3});
end
end
