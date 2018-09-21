%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EC_set = getECstring(EC_set,ecNumbers)
% 
% Provides a single string with all the EC numbers associated to a given
% protein with according to the next format: "ECX.X.X.X ECX.X.X.X ECX.X.X.X"
%
% Ivan Domenzain.      Last edited: 2018-09-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EC_set = getECstring(EC_set,ecNumbers)
%In case of several ec numbers for the same protein
new_EC_set = strsplit(ecNumbers,' ');
for l = 1:length(new_EC_set)
    EC_set = [EC_set 'EC' new_EC_set{l}];
    if l<length(new_EC_set)
        EC_set = [EC_set ' '];
    end
end
end
