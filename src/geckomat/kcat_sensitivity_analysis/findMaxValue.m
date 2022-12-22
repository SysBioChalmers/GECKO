%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [value,organism,parameter] = findMaxValue(EC_cell,BRENDA, SA_cell)
%
% Function that gets the maximum kinetic parameter (Kcat or S.A.*Mw) from 
% the BRENDA files for the specified set of EC numbers. The algorithm also 
% returns the organism and the parameter type (Kcat or S.A.) of the query.
%
% Ivan Domenzain    Last edited. 2018-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,organism,parameter] = findMaxValue(EC_cell,BRENDA,SA_cell)
    %Looks for the maximum turnover number available for the EC# associated
    %with the uniprot code
    EC_cell    = strsplit(EC_cell,' ');
    value      = [];
    organism   = [];
    parameter  = [];
    for i=1:length(EC_cell)
        find_flag  = false;
        %In case that wild cards are present in the EC number the search on
        %the BRENDA file will be relaxed.
        if ~isempty(strfind(EC_cell{i},'-'))
             EC_cell{i} = EC_cell{i}(strfind(EC_cell{i},'-')-1:end);
             find_flag  = true;
        end    
        ECnumber = ['EC' EC_cell{i}];
        Kcat     = 0; orgK = '';
        if find_flag == true
            matching = find(~cellfun(@isempty,strfind(BRENDA{1},ECnumber)));
        else
            % If no wild cards are present the EC number search in the
            % BRENDA file will look for an exact match
            matching = find(strcmpi(ECnumber,BRENDA{1}));
        end
        %Gets the maximum Kcat value for the queried EC#
        if ~isempty(matching)
            [Kcat, maxIndx] = max(BRENDA{4}(matching));
            orgK            = BRENDA{3}(matching(maxIndx));
        end        
        % Looks for the maximum SA*Mw value available for the EC number
        SA_Mw = 0; orgS = '';
        if find_flag == true
            matching = find(~cellfun(@isempty,strfind(SA_cell{1},ECnumber)));
        else
            matching = find(strcmpi(ECnumber,SA_cell{1}));
        end
        %Gets the maximum SA*Mw value for the queried EC#
        if ~isempty(matching)
            [SA_Mw, maxIndx] = max(SA_cell{3}(matching));
            SA_Mw            = SA_Mw; %[1/hr]
            orgS             = SA_cell{2}(matching(maxIndx));
        end        
        %Choose the maximal available value as a turnover number for the EC
        value  = [value; max(Kcat,SA_Mw)]; 

        if Kcat>SA_Mw
            organism  = [organism; {orgK}];
            parameter = [parameter; {'K_cat'}];
        else
        	organism  = [organism; {orgS}];
            parameter = [parameter; {'SA*Mw'}];
        end   
    end
    [value, maxIndx] = max(value);
    organism         = organism(maxIndx);
    parameter        = parameter(maxIndx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
