function model = findMetSmiles(model, modelAdapter, verbose)
% findMetSMILES
%   Queries PubChem by metabolite names to obtain SMILES. Matches will also
%   be stored in userData/***/data/smilesDB.tsv, that will also be queried
%   first next time the function is run. If the model already has a
%   metSmiles field, then non-empty entries will not be overwritten.
%
% Input:
%   model        Input model, whose model.metNames field is used to find the
%                relevant SMILES
%   modelAdapter a loaded model adapter (Optional, will otherwise use the
%                default model adapter).
%   verbose      logical whether progress should be reported (Optional,
%                default true)
% Ouput:
%   model       Output model with model.metSmiles specified.
%
if nargin < 3 || isempty(verbose)
    verbose = true;
end
if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

[uniqueNames, ~, uniqueIdx] = unique(model.metNames);
uniqueSmiles(1:numel(uniqueNames),1) = {''};
protMets = startsWith(uniqueNames,'prot_');
metMatch = false(length(uniqueNames),1);
if verbose; fprintf('Check for local SMILES database... '); end
smilesDBfile = (fullfile(params.path,'data','smilesDB.tsv'));
if exist(smilesDBfile,'file')==2
    fID = fopen(smilesDBfile,'r');
    raw = textscan(fID,'%s %s','Delimiter','\t','HeaderLines',0);
    fclose(fID);
    smilesDB.names = raw{1};
    smilesDB.smile = raw{2};
    [metMatch, metIdx] = ismember(uniqueNames,smilesDB.names);
    uniqueSmiles(metMatch) = smilesDB.smile(metIdx(metMatch));
    if verbose; fprintf('done.\n'); end
else
    if verbose; fprintf('not found.\n'); end
end

if any(~metMatch & ~protMets)
    if verbose; fprintf('Querying PubChem for SMILES by metabolite names...   0%% complete'); end
    numUnique = numel(uniqueNames);
    for i = 1:numel(uniqueNames)
        if metMatch(i) || protMets(i)
            break;
        end
        if verbose && rem(i-1,floor(numUnique/100+1)) == 0
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
    if verbose; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n'); end
end
newSmiles = uniqueSmiles(uniqueIdx);
if ~isfield(model,'metSmiles') || all(cellfun(@isempty,model.metSmiles))
    model.metSmiles = newSmiles;
else
    emptySmiles = cellfun(@isempty,model.metSmiles);
    model.metSmiles(emptySmiles) = newSmiles(emptySmiles);
end
out = [uniqueNames(~protMets), uniqueSmiles(~protMets)]';
fID = fopen(smilesDBfile,'w');
fprintf(fID,'%s\t%s\n',out{:});
fclose(fID);
end
