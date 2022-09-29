function metSmiles = findMetSmiles(metNames)
% findMetSMILES
%   Queries PubChem by metabolite names to obtain SMILES.
%
% Input:
%   metNames    array with metabolite names (e.g. model.metNames)
%
% Ouput:
%   metSmiles   array with matching SMILES. If model.metNames is used as
%               input, metSmiles can be specified as model.metSmiles.
% 

[uniqueNames, ~, uniqueIdx]          = unique(metNames);
uniqueSmiles(1:numel(uniqueNames),1) = {''};

fprintf('Querying PubChem for SMILES by metabolite names...   0%% complete');
numUnique = numel(uniqueNames);
for i = 1:numel(uniqueNames)
    if rem(i-1,floor(numUnique/100)) == 0
        progress = num2str(floor(100*(i/numUnique)));
        progress = pad(progress,3,'left');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
    end
    try
        uniqueSmiles{i,1} = strtrim(webread(['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' uniqueNames{i} '/property/CanonicalSMILES/TXT']));
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n');

metSmiles = uniqueSmiles(uniqueIdx);
end