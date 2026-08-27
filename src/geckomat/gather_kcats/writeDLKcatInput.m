function writtenTable = writeDLKcatInput(model, varargin)
% writeDLKcatInput  Prepare and write the input file for DLKcat.
%
% Prepares the input for DLKcat, and writes it to data/DLKcat.tsv in the
% obj.params.path specified in the ModelAdapter.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% ecRxns : logical
%     for which reactions (from model.ec.rxns) DLKcat should predict kcat
%     values, provided as logical vector with same length as model.ec.rxns
%     (default: all reactions).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% onlyWithSmiles : logical
%     whether to only include metabolites with SMILES (default true).
% filename : char
%     path to the input file, including the filename and .tsv extension
%     (default: data/DLKcat.tsv from the obj.params.path folder specified in
%     the modelAdapter).
% overwrite : logical
%     whether existing file should be overwritten (default false, to prevent
%     overwriting a file that already contains DLKcat-predicted kcat values).
%
% Returns
% -------
% writtenTable : cell
%     the table written, mainly to be used for testing purposes.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     writtenTable = writeDLKcatInput(model);
%     writtenTable = writeDLKcatInput(model, 'onlyWithSmiles', false);
%
% See also
% --------
% runDLKcat, readDLKcatOutput

[geckoPath, ~] = findGECKOroot();

p = parseGECKOargs(varargin, { ...
    'ecRxns',         []; ...
    'modelAdapter',   []; ...
    'onlyWithSmiles', []; ...
    'filename',       []; ...
    'overwrite',      []});
ecRxnsSel      = p.ecRxns;
modelAdapter   = p.modelAdapter;
onlyWithSmiles = p.onlyWithSmiles;
filename       = p.filename;
overwrite      = p.overwrite;

if isempty(ecRxnsSel)
    ecRxnsSel = true(numel(model.ec.rxns),1);
elseif ~logical(ecRxnsSel)
    error('ecRxns should be provided as logical vector')
elseif numel(ecRxnsSel)~=numel(model.ec.rxns)
    error('Length of ecRxns is not the same as model.ec.rxns')
end
ecRxnsSel = find(ecRxnsSel); % Change to indices, into the full model.ec.rxns

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(onlyWithSmiles)
    onlyWithSmiles=true;
end

if isempty(filename)
    filename = fullfile(params.path,'data','DLKcat.tsv');
elseif ~endsWith(filename,'.tsv')
    error('If filename is provided, it should include the .tsv extension.')
end

if isempty(overwrite) || ~overwrite % If is true
    if exist(filename,'file')
        error([filename ' already exists, either delete it first, or set the overwrite input argument as true'])
    end
end

if ~model.ec.geckoLight
   origRxns = model.ec.rxns;
else
   origRxns = extractAfter(model.ec.rxns,4);
end
origRxnsToInclude = origRxns(ecRxnsSel);

% Map back to original reactions, to extract substrates
[sanityCheck,origRxnIdxs] = ismember(origRxnsToInclude,model.rxns);
if ~all(sanityCheck)
    error('Not all reactions in model.ec.rxns are found in model.rxns')
end

% Ignore selected metabolites (metal ions, proteins etc.). First check by
% name (case insensitive, without white spaces and special characters),
% then also try to match with metSmiles (if available).
metsNoSpecialChars = lower(regexprep(model.metNames,'[^0-9a-zA-Z]+',''));
if exist(fullfile(params.path,'data','DLKcatIgnoreMets.tsv'),'file')
    fID        = fopen(fullfile(params.path,'data','DLKcatIgnoreMets.tsv'));
else
    fID        = fopen(fullfile(geckoPath,'databases','DLKcatIgnoreMets.tsv'));
end
fileData   = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[ignoreMets, ignoreSmiles] = deal(fileData{[1,2]});
ignoreMets = lower(regexprep(ignoreMets,'[^0-9a-zA-Z]+',''));
ignoreSmiles(cellfun(@isempty,ignoreSmiles)) = [];

ignoreMetsIdx  = logical(ismember(metsNoSpecialChars,ignoreMets));
if isfield(model,'metSmiles')
    ignoreMetsIdx = ignoreMetsIdx | logical(ismember(model.metSmiles,ignoreSmiles));
end
% Also leave out protein-usage pseudometabolites
ignoreMetsIdx = ignoreMetsIdx | startsWith(model.mets,'prot_');
reducedS = model.S;
reducedS(ignoreMetsIdx,:) = 0;

% Ignore currency metabolites if they occur in pairs. First check by
% name (case insensitive, without white spaces and special characters),
% then also try to match with metSmiles (if available).
if exist(fullfile(params.path,'data','DLKcatCurrencyMets.tsv'),'file')
    fID = fopen(fullfile(params.path,'data','DLKcatCurrencyMets.tsv'));
else
    fID = fopen(fullfile(geckoPath,'databases','DLKcatCurrencyMets.tsv'));
end
fileData = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[currencyMets(:,1), currencyMets(:,2)] = deal(fileData{[1,2]});
currencyMets = lower(regexprep(currencyMets,'[^0-9a-zA-Z]+',''));

for i=1:size(currencyMets,1)
    subs = strcmp(currencyMets(i,1),metsNoSpecialChars);
    prod = strcmp(currencyMets(i,2),metsNoSpecialChars);
    [~,subsRxns]=find(reducedS(subs,:));
    [~,prodRxns]=find(reducedS(prod,:));
    pairRxns = intersect(subsRxns,prodRxns);
    tempRedS=reducedS;
    tempRedS([find(subs);find(prod)],pairRxns) = 0;
    % Do not remove currency mets if no substrate remains
    rxnsWithRemainingSubstrates = any(tempRedS(:,pairRxns) < 0,1);
    reducedS([find(subs);find(prod)],intersect(pairRxns,pairRxns(rxnsWithRemainingSubstrates))) = 0;
end

% origRxnIdxs already selects exactly the model.rxns columns for the requested
% ecRxns subset (one column per requested reaction, via origRxns(ecRxns) above), so
% clearedRedS needs no further filtering here.
clearedRedS = reducedS(:,origRxnIdxs);

% Enumerate all substrates for each reaction. clearedRedS has one column per entry
% of ecRxnsSel, in that same order, so `reactions` here indexes *that subset*, not
% model.ec.rxns directly --- translate back to a real model.ec.rxns index before
% using it against model.ec.rxns / model.ec.rxnEnzMat, which are always full-length.
[substrates, reactions] = find(clearedRedS<0);
globalEcRxns = ecRxnsSel(reactions);

% Enumerate all proteins for each reaction
[proteins, enzRows] = find(transpose(model.ec.rxnEnzMat(globalEcRxns,:)));

% Prepare output
out(1,:) = model.ec.rxns(globalEcRxns(enzRows));
out(2,:) = model.ec.genes(proteins);
out(3,:) = model.metNames(substrates(enzRows));
if isfield(model,'metSmiles')
    out(4,:) = model.metSmiles(substrates(enzRows));
else
    out(4,:) = cell(numel(substrates(enzRows)),1);
end

out(5,:) = model.ec.sequence(proteins);
if onlyWithSmiles
    out(:,cellfun(@isempty,out(4,:))) = [];
else
    out(4,cellfun(@isempty,out(4,:))) = {'None'};
end
out(6,:) = cell(numel(out(1,:)),1);
out(6,:) = {'NA'};

% Write file
fID = fopen(filename,'w');
fprintf(fID,'%s\t%s\t%s\t%s\t%s\t%s\n',out{:});
fclose(fID);
fprintf('Model-specific DLKcat input stored at %s\n',filename);

writtenTable = out;
end
