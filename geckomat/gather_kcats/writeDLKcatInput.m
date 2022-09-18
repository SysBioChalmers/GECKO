function DLKcatIDs = writeDLKcatInput(model,inFile,ecRxns)
% writeDLKcatInput
%   Prepares the input file that can be used by DLKcat
%
% Input:
%   model       an ec-model in RAVEN format
%   inFile      name and path of the DLKcat input file to be written.
%               (Opt, default is 'DLKcatInput.tsv' in the current working
%               directory)
%   ecRxns      for which reactions (from model.ec.rxns) DLKcat should
%               predict kcat values, provided as logical vector with same
%               length as model.ec.rxns. (Opt, default is all reactions)
%
% Output:
%   DLKcatIDs   vector that specify reaction identifiers and genes and for
%               each line in the DLKcat input file that is written. To be
%               used for matching the DLKcat results back to the model.
%

if nargin<2
    inFile='DLKcatInput.tsv';
elseif endsWith(inFile,filesep()) % Append file name if only path is given
    inFile=fullfile(inFile,'DLKcatInput.tsv');
end
if nargin<3
    ecRxns=true(numel(model.ec.rxns),1);
elseif ~logical(ecRxns)
    error('ecRxns should be provided as logical vector')
elseif numel(ecRxns)~=numel(model.ec.rxns)
    error('Length of ecRxns is not the same as model.ec.rxns')
end
[geckoPath, ~] = findGECKOroot();

% Identify reactions for which kcat should be predicted (entry in model.ec.rxns)
rxnsToInclude=model.ec.rxns(ecRxns);
ecRxns=find(ecRxns); % Change to indices
% Map back to original reactions, to extract substrates
[sanityCheck,rxnIdxs]=ismember(rxnsToInclude,model.rxns);
if ~all(sanityCheck)
    error('Not all reactions in model.ec.rxns are found in model.rxns')
end

% Specify currency metabolites that should be ignored as substrates
%TODO: rename to ignoreMets, can also include e.g. large metabolites
fID          = fopen(fullfile(geckoPath,'databases','smallMets.tsv'));
smallMets    = textscan(fID,'%s','delimiter','\t');
fclose(fID);
fID          = fopen(fullfile(geckoPath,'databases','currencyMets.tsv'));
currencyMets = textscan(fID,'%s','delimiter','\t');
fclose(fID);
%TODO: currencyMets check as substrate and product
ignoreMets   = [smallMets{1}; currencyMets{1}];
ignoreMets   = find(ismember(model.metNames,ignoreMets));
% Enumerate all substrates for each reaction
[substrates, reactions] = find(model.S(:,rxnIdxs)<0);
toIgnore                = ismember(substrates,ignoreMets);
substrates(toIgnore)    = [];
reactions(toIgnore)     = [];
% Enumerate all proteins for each reaction
[proteins,ecRxns] = find(transpose(model.ec.rxnEnzMat(ecRxns(reactions),:)));

% Prepare output
clear out
out(1,:) = model.metNames(substrates(ecRxns));
if isfield(model,'metSmiles')
    out(2,:) = model.metSmiles(substrates(ecRxns));
else
    out(2,:) = cell(numel(substrates(ecRxns)),1);
end
out(2,cellfun(@isempty,out(2,:))) = {'None'};
out(3,:) = model.ec.sequence(proteins);
out = [{'Substrate Name';'Substrate SMILES';'Protein Sequence'}, out];

% Write file
fID = fopen(inFile,'w');
fprintf(fID,'%s\t%s\t%s\n',out{:});
fclose(fID);

% Define DLKcatIDs output
DLKcatIDs.rxns       = model.rxns(rxnIdxs(reactions(ecRxns)));
DLKcatIDs.genes      = model.ec.genes(proteins);
DLKcatIDs.substrates = transpose(out(1,2:end));
end
