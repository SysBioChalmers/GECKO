%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [KCATcell, SAcell] = loadBRENDAdata 
     % TODO: potential refactor with DLKcat derived BRENDA data
     geckoPath      = findGECKOroot;
     KCAT_file      = fullfile(geckoPath,'databases','max_KCAT.txt');
     SA_file        = fullfile(geckoPath,'databases','max_SA.txt');
     MW_file        = fullfile(geckoPath,'databases','max_MW.txt');
    
     %Extract BRENDA DATA from files information
     scalingFactor = 1/60;    %[umol/min/mg] -> [mmol/s/g]    Old: 60 [umol/min/mg] -> [mmol/h/g]
     SA            = openDataFile(SA_file,scalingFactor);
     scalingFactor = 1/1000;  %[g/mol] -> [g/mmol]
     MW            = openDataFile(MW_file,scalingFactor);
     
     for i=1:4
         SAcell{i} = [];
     end
     previousEC = []; EC_indexes = [];
     
    %build an index on MW{1} to speed things up a bit
    %first just extract the genus (i.e. the first part of the name)
    MWECNum = upper(unique(MW{1}));
    MWECNumIndices = cell(length(MWECNum),1);
    MWECNumHashMap = containers.Map(MWECNum,1:length(MWECNum));
    for i = 1:length(MW{1})
         matchInd = cell2mat(values(MWECNumHashMap, MW{1}(i)));
         MWECNumIndices{matchInd} = [MWECNumIndices{matchInd};i];
    end

     
     for i=1:length(SA{1})
         %Gets the indexes of the EC repetitions in the MW cell for every
         %new (different) EC
         if ~strcmpi(SA{1}(i), previousEC)
             key = upper(SA{1}(i));
             if isKey(MWECNumHashMap,key) %annoyingly, this seems to be needed
                 matchInd = cell2mat(values(MWECNumHashMap,key));
                 EC_indexes = MWECNumIndices{matchInd};
             else
                 EC_indexes = [];
             end
         end
         mwEC{1} = MW{3}(EC_indexes); mwEC{2} = MW{4}(EC_indexes);
         % just looks for the first match because just the maximal value for
         % each EC# / Orgaism is reported on the file
         org_index = find(strcmpi(SA{3}(i),mwEC{1}),1);
         if ~isempty(org_index)
             SAcell{1} = [SAcell{1};SA{1}(i)];
             SAcell{2} = [SAcell{2};SA{3}(i)];
             SAcell{3} = [SAcell{3}; SA{4}(i)*mwEC{2}(org_index)]; %[1/hr]
             SAcell{4} = [SAcell{4}; mwEC{2}(org_index)];
         end
         previousEC = SA{1}(i);
     end
     
     %remove EC in front of all the EC numbers
     if ~isempty(KCATcell{1})
         KCATcell{1} = extractAfter(KCATcell{1},2);
     end
     if ~isempty(SAcell{1})
        SAcell{1} = extractAfter(SAcell{1},2);
     end

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function string_cells = stringSplit(cell_array)
         string_cells = {textscan(cell_array,'%q','Delimiter','//')};
         string_cells = string_cells{1}(1);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function data_cell = openDataFile(fileName,scalingFactor)
     fID          = fopen(fileName);
     data_cell    = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
     fclose(fID);
     data_cell{4} = data_cell{4}*scalingFactor;
     %Split string for each organism in the BRENDA data 
     %{name, taxonomy, KEGG code}
     data_cell{3}  = cellfun(@stringSplit, data_cell{3});
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
