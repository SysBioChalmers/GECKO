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
    if rem(i-1,floor(numUnique/100+1)) == 0
        progress = num2str(floor(100*(i/numUnique)));
        progress = pad(progress,3,'left');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
    end
        retry = true;
        while retry
            try
                retry = false;
                smileResult       = webread(['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' uniqueNames{i} '/property/CanonicalSMILES/TXT']);
                %Sometimes multiple lines are given, with alternative SMILES. Only
                %keep the first suggestion.
                smileResult       = regexp(smileResult,'(^\S*)\n','once','tokens');
                uniqueSmiles{i,1} = smileResult{1,1};
            catch exception
                %Sometimes the call fails, for example since the server is busy. In those cases
                %we should retry until we get a response. Some errors however are because the metabolite
                %name doesn't exist in the database (404) or some other error (the metabolite contains
                %a slash or similar, 400) - in those cases we need to give up, otherwise the function
                %will enter an infinite loop.
                if ~(strcmp(exception.identifier, 'MATLAB:webservices:HTTP404StatusCodeError') || ...
                     strcmp(exception.identifier, 'MATLAB:webservices:HTTP400StatusCodeError'))
                    retry = true;
                end
            end
        end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n');

metSmiles = uniqueSmiles(uniqueIdx);
end