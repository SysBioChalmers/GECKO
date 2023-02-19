function writtenTable = writeDLKcatInput(model, ecRxns, modelAdapter, onlyWithSmiles, filename)
% writeDLKcatInput
%   Prepares the input for DLKcat, and writes it to data/DLKcat.tsv
%   in the obj.params.path specified in the ModelAdapter.
%
% Input:
%   model           an ec-model
%   ecRxns          for which reactions (from model.ec.rxns) DLKcat should
%                   predict kcat values, provided as logical vector with
%                   same length as model.ec.rxns. (Opt, default is all
%                   reactions)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%   onlyWithSmiles  logical whether to only include metabolites with SMILES
%                   (optional, default true)
%   filename        Filename (Optional). Normally this parameter should not be 
%                   supplied, but it is useful for test cases.
%
% Output:
%   writtenTable    The table written, mainly to be used for testing purposes.

[geckoPath, ~] = findGECKOroot();

if nargin<2 || isempty(ecRxns)
    ecRxns = true(numel(model.ec.rxns),1);
elseif ~logical(ecRxns)
    error('ecRxns should be provided as logical vector')
elseif numel(ecRxns)~=numel(model.ec.rxns)
    error('Length of ecRxns is not the same as model.ec.rxns')
end
ecRxns = find(ecRxns); % Change to indices

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin<4 || isempty(onlyWithSmiles)
    onlyWithSmiles=true;
end

if nargin<5 || isempty(filename)
    filename = fullfile(params.path,'data','DLKcat.tsv');
end

if ~model.ec.geckoLight
   origRxns = model.ec.rxns;
else
   origRxns = extractAfter(model.ec.rxns,4);
end
origRxnsToInclude = origRxns(ecRxns);

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

ignoreMetsIdx  = logical(ismember(metsNoSpecialChars,ignoreMets));
if isfield(model,'metSmiles')
    ignoreMetsIdx = ignoreMetsIdx | logical(ismember(model.metSmiles,ignoreSmiles));
end
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
fileData = textscan(fID,'%s %s %s %s','delimiter','\t');
fclose(fID);
[currencyMets(:,1), currencyMets(:,2), currSmiles(:,1), currSmiles(:,2)] = deal(fileData{[1,3,2,4]});
[currencyMets(:,1), currencyMets(:,2)] = deal(fileData{[1,2]});
currencyMets = lower(regexprep(currencyMets,'[^0-9a-zA-Z]+',''));

for i=1:size(currencyMets,1)
    subs = strcmp(currencyMets(i,1),metsNoSpecialChars);
    prod = strcmp(currencyMets(i,2),metsNoSpecialChars);
    if isfield(model,'metSmiles')
        subs = subs | logical(ismember(model.metSmiles,currSmiles(i,1)));
        prod = prod | logical(ismember(model.metSmiles,currSmiles(i,2)));
    end
    [~,subsRxns]=find(reducedS(subs,:));
    [~,prodRxns]=find(reducedS(prod,:));
    pairRxns = intersect(subsRxns,prodRxns);
    reducedS([subs;prod],pairRxns) = 0;
end

%filter out the reactions we're not interested in - will solve the problem for both full and light
clearedRedS = reducedS(:,origRxnIdxs);
rxnsToClear = true(length(origRxnIdxs),1);
rxnsToClear(ecRxns) = false;
clearedRedS(:,rxnsToClear) = 0;

% Enumerate all substrates for each reaction
[substrates, reactions] = find(clearedRedS<0); %the reactions here are in model.ec.rxns space

% Enumerate all proteins for each reaction
[proteins, ecRxns] = find(transpose(model.ec.rxnEnzMat(reactions,:)));

% Prepare output
out(1,:) = model.ec.rxns(reactions(ecRxns));
out(2,:) = model.ec.genes(proteins);
out(3,:) = model.metNames(substrates(ecRxns));
if isfield(model,'metSmiles')
    out(4,:) = model.metSmiles(substrates(ecRxns));
else
    out(4,:) = cell(numel(substrates(ecRxns)),1);
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

writtenTable = out;
end
