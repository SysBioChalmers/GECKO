%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ecModel = flexibilizeKcats(ecModel,ecModel_batch,PDB,limKcats,gRexp)
%
% Function that gets the limiting Kcat values in an EC model (according to
% a sensitivity analysis), then it modifies each of those values according to 
% the maximal available values in the BRENDA files (Kcats and SA*Mw) when a
% manual curated option is not specified.
% 
% The algorithm iterates until the model grows at the same rate provided 
% by the user (batch growth on glucose minimal media recommended)
%
% Ivan Domenzain    Last edited. 2018-02-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ecModel = flexibilizeKcats(ecModel,ecModelBatch,gRexp)
    

    modified_kcats = []; modifiedRxns = [];
    modifications  = []; error        = [];  

    e  = -100; i=1; 
    current = pwd; 
    %Load BRENDA data:
    KCAT_file          = 'max_KCAT.txt';
    SA_file            = 'max_SA.txt';
    MW_file            = 'max_MW.txt';
    [BRENDA,SA_cell]   = loadBRENDAdata(KCAT_file,SA_file,MW_file);
    
    %Iterates while growth rate is being underpredicted
    disp('********************Limiting Kcats curation********************')
    while e<=0
        cd (current)
        %Get the top growth rate-limiting enzyme (uniprot code basis)
        [limKcat,limRxns,breakFlag] =findTopLimitations(ecModelBatch,modified_kcats);
        
        if breakFlag == false
            [ecModelBatch,data] = changeKcat(ecModelBatch,limKcat,gRexp,...
                                                         BRENDA,SA_cell);
                                          
            cd (current)
            %Saves the parameter modification information
            modifications = [modifications; data];
            e             = data{1,7};
            %Add a string with the uniprot code and the rxn number in order
            %to keep track of the modified coefficients
            str            = {horzcat(data{1},'_',num2str(limKcat{3}))};
            modified_kcats = [modified_kcats; str];
            disp([data{1} '_' num2str(limKcat{3}) ': ' limKcat{6}(1)])
            
            error  = [error; e];           
            str    = ['#' num2str(i) ' prev:' num2str(data{1,7}) ' new:' ...
                      num2str(data{1,8}) ' CC:' num2str(limKcat{1,5}) ...
                                                       ' Err:' num2str(e)];
            disp(str)          
            i = i+1;            
        else
            break
        end
    end  
    cd (current)
    if ~isempty(modifications)
        [m,n]         = size(ecModel.S);
        ecModel.S     = ecModelBatch.S(1:m,1:n);
        varNamesTable = {'Unicode','enz_pos','rxn_pos','Organism',...
                         'Modified','Parameter','oldValue','newValue',...
                         'error','gRControlCoeff'}:

        modifications = cell2table(modifications,'VariableNames',varNamesTable);
        writetable(modifications, 'KcatModifications.txt');
    else
        if ~isempty(limRxns)
            varNamesTable = {'rxnNames','rxnPos','gRControlCoeff'};
            modifications = cell2table(modifications,'VariableNames',varNamesTable);
            writetable(modifications, 'LimitingRxns.txt');
        end
    end
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,data] = changeKcat(model,limKcats,gR_exp,BRENDA,SA_cell)
   % Gets the Unicode 
    UniCode  = limKcats{1}{1}(strfind(limKcats{1}{1},'_')+1:end);
    % Map the UNIPROT code (kcat)
    [ECnumber, Mw] = findECnumber(UniCode);
    enzIndx = limKcats{2}(1);
    rxnIndx = limKcats{3}(1);
    data    = {UniCode,[],[],[],[],[],[],[]};
    e       = 0;
    
    if ~isempty(ECnumber) && ~isempty(Mw)
    	flag           = false;
    	previous_value = -1/(3600*model.S(enzIndx,rxnIndx)); %[1/s]
        
            disp(['Automatic search // ' 'EC#: ' ECnumber])
           %Looks for the maximal value available for the respective EC
           %number (Kcats and SA*Mw if available)
            [Kcat,org, match] = findMaxValue(ECnumber,BRENDA,SA_cell,Mw);
            coeff             = -1/(Kcat);  
           %Change the kinetic coefficient just if a higher value was found
            if coeff > model.S(enzIndx,rxnIndx)
            	flag = true;
                model.S(enzIndx,rxnIndx) = coeff;
            end
        new_value = -1/(3600*model.S(enzIndx,rxnIndx));
           
        % After changing the i-th kcat limiting value a simulation is
        % performed and the growth rate and absolute error are saved 
        model_sim            = model;
        gR_pos               = find(strcmpi(model_sim.rxnNames,'growth'));
        model_sim.c          = zeros(size(model_sim.c));
        model_sim.c(gR_pos)  = 1;
        solution             = solveLP(model_sim);
        model_sim.lb(gR_pos) = 0.999*solution.x(gR_pos);
        model_sim.ub(gR_pos) = solution.x(gR_pos);
        solution             = solveLP(model_sim,1);
        e                    = ((solution.x(gR_pos)-gR_exp)/gR_exp)*100;
        
        data                 = {UniCode,limKcats{2}(1),limKcats{3}(1),...
                                org,flag,match,previous_value,new_value,...
                                                         e,limKcats{5}(1)};
   end 
         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [ECnumber, Mw] = findECnumber(Unicode)
    cd ../../Databases
    load ('ProtDatabase.mat')
    DB1{1} = swissprot(:,1);DB1{2} = swissprot(:,4);DB1{3} = swissprot(:,5);
    DB2{1} = kegg(:,1);     DB2{2} = kegg(:,4);     DB2{3} = kegg(:,5);
    ECnumber = {};
    % First look for the UNIPROT ID in the swissprot DB structure
    matching = indexes_string(DB1{1},Unicode,true);  
    if ~isempty(matching)
        ECnumber = DB1{2}{matching};
        Mw       = DB1{3}{matching};
    end
    % If nothing comes up then look into the KEGG DB structure
    if isempty(ECnumber)
        matching = indexes_string(DB2{1},Unicode,true);
        if ~isempty(matching)
            ECnumber = DB2{2}{matching};
            Mw       = DB1{3}{matching};
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matching = indexes_string(cell_array,str,flag)
    matching  = strfind(cell_array,str);
    if flag == false
        matching = find(~cellfun(@isempty,matching));
    else
        matching = find(~cellfun(@isempty,matching),1);
    end
end
