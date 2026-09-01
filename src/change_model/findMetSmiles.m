function [model,noSMILES] = findMetSmiles(model, varargin)
% findMetSmiles  Query PubChem by metabolite names to obtain SMILES.
%
% Matches are cached in tutorials/***/data/smilesDB.tsv, which is checked
% first on subsequent runs so repeated queries are avoided. If the model
% already has a metSmiles field, non-empty entries are left unchanged.
%
% Parameters
% ----------
% model : struct
%     a model whose metNames field is used to find the relevant SMILES.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% verbose : logical
%     whether progress should be reported (default true).
%
% Returns
% -------
% model : struct
%     model with model.metSmiles specified.
% noSMILES : cell
%     metabolite names for which no SMILES could be found.
%
p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'verbose',      true});
modelAdapter = p.modelAdapter;
verbose      = p.verbose;
if isempty(verbose); verbose = true; end

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

[uniqueNames, ~, uniqueIdx] = unique(regexprep(model.metNames,'^prot_.*',''));
uniqueSmiles(1:numel(uniqueNames),1) = {''};
metMatch = false(length(uniqueNames),1);
metMatch(strcmp(uniqueNames,'')) = 1; % No need trying to match empty fields
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

if any(~metMatch)
    PB = progressReport(numel(uniqueNames), 'Querying PubChem for SMILES by metabolite names');
    webOptions = weboptions('Timeout', 30);
    for i = 1:numel(uniqueNames)
        PB.count;
        if metMatch(i)
            continue;
        end
        retry = 0;
        while retry < 10
            try
                smileResult       = webread(['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' uniqueNames{i} '/property/CanonicalSMILES/TXT'], webOptions);
                %Sometimes multiple lines are given, with alternative SMILES. Only
                %keep the first suggestion.
                smileResult       = regexp(smileResult,'(^\S*)\n','once','tokens');
                uniqueSmiles{i,1} = smileResult{1,1};
                retry = 15; % success: no retry
            catch exception
                %Sometimes the call fails, for example since the server is busy. In those cases
                %we will try 10 times. Some errors however are because the metabolite
                %name does no exist in the database (404) or some other error (the metabolite contains
                %a slash or similar, 400 or 500). In those cases we can
                %immediately give up.
                if (strcmp(exception.identifier, 'MATLAB:webservices:HTTP404StatusCodeError') || ...
                        strcmp(exception.identifier, 'MATLAB:webservices:HTTP400StatusCodeError') || ...
                        strcmp(exception.identifier, 'MATLAB:webservices:HTTP500StatusCodeError'))
                    uniqueSmiles(i) = {''};
                    retry = 15;
                else
                    retry = retry + 1;
                end
            end
        if retry == 10
            error('Cannot reach PubChem. Check your internet connection and try again.')
        end
        end
        % Append one line each time, in case internet connection is lost
        % halfway. Open & close file each time to avoid leaving the file
        % open when breaking the function.
        out = [uniqueNames(i), uniqueSmiles(i)];
        fID = fopen(smilesDBfile,'a');
        fprintf(fID,'%s\t%s\n',out{:});
        fclose(fID);
    end
    PB.done;
    if verbose
        fprintf('Model-specific SMILES database stored at %s\n',smilesDBfile);
    end
end
newSmiles = uniqueSmiles(uniqueIdx);
noSMILES = cellfun(@isempty,uniqueSmiles);
successRatio = 1-(numel(find(noSMILES))/numel(uniqueSmiles));
fprintf('SMILES could be found for %s%% of the unique metabolite names.\n',num2str(successRatio*100,'%.0f'))
noSMILES = uniqueNames(noSMILES);

if ~isfield(model,'metSmiles') || all(cellfun(@isempty,model.metSmiles))
    model.metSmiles = newSmiles;
else
    emptySmiles = cellfun(@isempty,model.metSmiles);
    model.metSmiles(emptySmiles) = newSmiles(emptySmiles);
end
end
