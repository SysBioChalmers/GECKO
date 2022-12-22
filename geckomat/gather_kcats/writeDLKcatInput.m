function writeDLKcatInput(model,inFile,ecRxns,ignoreMets,currencyMets)
% writeDLKcatInput
%   Prepares the input file that can be used by DLKcat
%
% Input:
%   model           an ec-model in RAVEN format
%   inFile          name and path of the DLKcat input file to be written.
%                   (Opt, default is 'DLKcatInput.tsv' in the current
%                   working directory)
%   ecRxns          logical of length model.ec.rxns that specifies for
%                   which reactions DLKcat should predict kcat values
%                   (optional, by default all model.ec.eccodes entries
%                   are populated by this function)
%   ignoreMets      cell array with metabolite names that should not be
%                   included in the DLKcat input file. (optional, default
%                   it loads GECKO/databases/DLKcatIgnoreMets.tsv)
%   currencyMets    cell array (with 2 columns) with pairs of currency 
%                   metabolites (identified by their metabolite name) that
%                   that should not be included in the DLKcat input file
%                   for reactions where both are involved. (optional,
%                   default it loads GECKO/databases/DLKcatIgnoreMets.tsv)
%

if nargin<2
    inFile = 'DLKcatInput.tsv';
elseif endsWith(inFile,filesep()) % Append file name if only path is given
    inFile = fullfile(inFile,'DLKcatInput.tsv');
end
if nargin<3
    ecRxns = true(numel(model.ec.rxns),1);
elseif ~logical(ecRxns)
    error('ecRxns should be provided as logical vector')
elseif numel(ecRxns)~=numel(model.ec.rxns)
    error('Length of ecRxns is not the same as model.ec.rxns')
end
[geckoPath, ~] = findGECKOroot();

% Identify reactions for which kcat should be predicted (entry in model.ec.rxns)
rxnsToInclude = model.ec.rxns(ecRxns);
ecRxns        = find(ecRxns); % Change to indices
% Map back to original reactions, to extract substrates
[sanityCheck,rxnIdxs] = ismember(rxnsToInclude,model.rxns);
if ~all(sanityCheck)
    error('Not all reactions in model.ec.rxns are found in model.rxns')
end

% Ignore selected metabolites (metal ions, proteins etc.). First check by
% name (case insensitive, without white spaces and special characters),
% then also try to match with metSmiles (if available).
metsNoSpecialChars = lower(regexprep(model.metNames,'[^0-9a-zA-Z]+',''));
if nargin<4
    fID        = fopen(fullfile(geckoPath,'databases','DLKcatIgnoreMets.tsv'));
    fileData   = textscan(fID,'%s %s','delimiter','\t');
    fclose(fID);
    [ignoreMets, ignoreSmiles] = deal(fileData{[1,2]});
    ignoreMets = lower(regexprep(ignoreMets,'[^0-9a-zA-Z]+',''));
end
ignoreMetsIdx  = logical(ismember(metsNoSpecialChars,ignoreMets));
if isfield(model,'metSmiles')
    ignoreMetsIdx = ignoreMetsIdx | logical(ismember(model.metSmiles,ignoreSmiles));
end
reducedS = model.S;
reducedS(ignoreMetsIdx,:) = 0;

% Ignore currency metabolites if they occur in pairs. First check by
% name (case insensitive, without white spaces and special characters),
% then also try to match with metSmiles (if available).
if nargin<5
    fID      = fopen(fullfile(geckoPath,'databases','DLKcatCurrencyMets.tsv'));
    fileData = textscan(fID,'%s %s %s %s','delimiter','\t');
    fclose(fID);
    [currMets(:,1), currMets(:,2), currSmiles(:,1), currSmiles(:,2)] = deal(fileData{[1,3,2,4]});
    [currMets(:,1), currMets(:,2)] = deal(fileData{[1,2]});
    currMets = lower(regexprep(currMets,'[^0-9a-zA-Z]+',''));
end
for i=1:size(currMets,1)
    subs = find(strcmp(currMets(i,1),metsNoSpecialChars));
    prod = find(strcmp(currMets(i,2),metsNoSpecialChars));
    [~,subsRxns]=find(reducedS(subs,:));
    [~,prodRxns]=find(reducedS(prod,:));
    pairRxns = intersect(subsRxns,prodRxns);
    reducedS([subs;prod],pairRxns) = 0;
end

% Enumerate all substrates for each reaction
[substrates, reactions] = find(reducedS(:,rxnIdxs)<0);
% Enumerate all proteins for each reaction
[proteins, ecRxns] = find(transpose(model.ec.rxnEnzMat(ecRxns(reactions),:)));

% Prepare output
out(1,:) = model.rxns(rxnIdxs(reactions(ecRxns)));
out(2,:) = model.ec.genes(proteins);
out(3,:) = model.metNames(substrates(ecRxns));
if isfield(model,'metSmiles')
    out(4,:) = model.metSmiles(substrates(ecRxns));
else
    out(4,:) = cell(numel(substrates(ecRxns)),1);
end
out(4,cellfun(@isempty,out(4,:))) = {'None'};
out(5,:) = model.ec.sequence(proteins);
out = [{'Reaction ID';'Gene ID';'Substrate Name';'Substrate SMILES';'Protein Sequence'}, out];

% Write file
fID = fopen(inFile,'w');
fprintf(fID,'%s\t%s\t%s\t%s\t%s\n',out{:});
fclose(fID);
end
