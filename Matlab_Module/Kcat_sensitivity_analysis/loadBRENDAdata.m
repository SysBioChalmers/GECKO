%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [KCATcell, SAcell] = loadBRENDAdata(KCAT_file,SA_file,MW_file )
     current = pwd;
     cd ../../Databases
     %Extract BRENDA DATA from files information
     scallingFactor = 3600;   %[1/s] -> [1/h]
     KCATcell      = openDataFile(KCAT_file,scallingFactor); 
     scallingFactor = 60;     %[umol/min/mg] -> [mmol/h/g]
     SA             = openDataFile(SA_file,scallingFactor); 
     scallingFactor = 1/1000; %[g/mol] -> [g/mmol]
     MW             = openDataFile(MW_file,scallingFactor); 
     
     for i=1:4
         SAcell{i} = [];
     end
     previousEC = []; EC_indexes = [];
     for i=1:length(SA{1})
         %Gets the indexes of the EC repetitions in the MW cell for every
         %new (different) EC
         if ~strcmpi(SA{1}(i), previousEC)
             EC_indexes = find(strcmpi(SA{1}(i),MW{1}));
         end
         mwEC{1} = MW{3}(EC_indexes); mwEC{2} = MW{4}(EC_indexes);
         % just looks for the first match because just the maximal value for
         % each EC# / Orgaism is reported on the file
         org_index = find(strcmpi(SA{3}(i),mwEC{1}),1);
         if ~isempty(org_index)
             SAcell{1} = [SAcell{1};SA{1}(i)];
             SAcell{2} = [SAcell{2};SA{3}(i)];
             SAcell{3} = [SAcell{3}; SA{4}(i)* mwEC{2}(org_index)]; %[1/hr]
             SAcell{4} = [SAcell{4}; mwEC{2}(org_index)];
         end
         previousEC = SA{1}(i);
     end
     cd (current)
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function data_cell = openDataFile(fileName,scallingFactor)
     fID          = fopen(fileName);
     data_cell    = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
     fclose(fID);
     data_cell{4} = data_cell{4}*scallingFactor;
     %Split string for each organism in the BRENDA data 
     %{name, taxonomy, KEGG code}
     data_cell{3}  = cellfun(@stringSplit, data_cell{3});
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function string_cells = stringSplit(cell_array)
        string_cells = {strsplit(cell_array,'//')};
        string_cells = string_cells{1}(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
