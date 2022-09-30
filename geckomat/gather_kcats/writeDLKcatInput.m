function writeDLKcatInput(model,inFile,ecRxns)
% writeDLKcatInput
%   Prepares the input file that can be used by DLKcat
%
% Input:
%   model           an ec-model in RAVEN format
%   inFile          name and path of the DLKcat input file to be written.
%                   (Opt, default is 'DLKcatInput.tsv' in the current
%                   working directory)
%   ecRxns          for which reactions (from model.ec.rxns) DLKcat should
%                   predict kcat values, provided as logical vector with
%                   same length as model.ec.rxns. (Opt, default is all
%                   reactions)
%   ignoreMets      cell array with metabolite names that should not be
%                   included in the DLKcat input file. (Opt, default it
%                   loads the list in GECKO/databases/DLKcatIgnoreMets.tsv)
%   currencyMets    cell array (with 2 columns) with pairs of currency 
%                   metabolites (identified by their metabolite name) that
%                   that should not be included in the DLKcat input file
%                   for reactions where both . (Opt, default it
%                   loads the list in GECKO/databases/DLKcatIgnoreMets.tsv)
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

% Ignore selected metabolites (metal ions, proteins etc.)
metsNoSpecialChars = lower(regexprep(model.metNames,'[^0-9a-zA-Z]+',''));
if nargin<4
    fID          = fopen(fullfile(geckoPath,'databases','DLKcatIgnoreMets.tsv'));
    ignoreMets   = textscan(fID,'%s','delimiter','\t');
    fclose(fID);
    ignoreMets   = lower(regexprep(ignoreMets{1},'[^0-9a-zA-Z]+',''));
end
ignoreMets   = logical(ismember(metsNoSpecialChars,ignoreMets));
reducedS     = model.S;
reducedS(ignoreMets,:) = 0;
% Ignore currency metabolites if they occur in pairs
if nargin<5
    fID          = fopen(fullfile(geckoPath,'databases','DLKcatCurrencyMets.tsv'));
    currencyMets = textscan(fID,'%s %s','delimiter','\t');
    fclose(fID);
    currencyMets    = lower(regexprep([currencyMets{1}, currencyMets{2}],'[^0-9a-zA-Z]+',''));
end
for i=1:size(currencyMets,1)
    subs = find(strcmp(currencyMets(i,1),metsNoSpecialChars));
    prod = find(strcmp(currencyMets(i,2),metsNoSpecialChars));
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
