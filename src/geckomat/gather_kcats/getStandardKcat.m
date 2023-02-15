function [model, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(model, modelAdapter, threshold)
% getStandardKcat
%   Calculate an standard kcat to rxns which a GPR have not been found in the
%   model
%
% Input:
%   model           an ecModel in GECKO 3 version
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%   threshold       a threshold to determine when use a kcat value based on
%                   the mean kcat of the reactions in the same subSystem or
%                   based on the median value of all the kcat in the model.
%                   Second option is used when the number of reactions in a
%                   determined subSystem is < threshold. (Optional, default = 10)
%
% Output:
%   model           ecModel where model.ec is expanded with a standard
%                   protein with standard kcat and standard MW, assigned to
%                   reactions without gene associations.
%   rxnsMissingGPR  a list of update rxns index with a standard value
%   standardMW      the standard MW value calculated
%   standardKcat    the standard Kcat value calculated
%
%   A pseudometabolite named 'prot_standard' will be added as well as a gene
%   and enzyme named 'standard'. While model.ec.kcat is populated,
%   applyKcatConstraints would still need to be run to apply the new
%   constraints to the S-matrix.
%
% Usage:
%    [model, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(model, modelAdapter, threshold);

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if nargin < 3
    threshold = 10;
end

% Maybe this can be an input ???
databases = loadDatabases('uniprot', modelAdapter);

% An stardard MW is defined for all the rxns which does not have a GPR
% rule defined. This is based in all the proteins reported for the specific
% organism in uniprot
standardMW = median(databases.uniprot.MW, 'omitnan');

% An standard Kcat is defined for all the rxns which does not have a GPR
% rule defined. In this case, the kcat value for a particular reaction is
% defined as the mean of the kcat values of the reactions involved in the
% same subsystem in which the given reaction is involved. Nevertheless, if
% a subSystem have a number of reactions lower than a treshold, the kcat
% value will be the median of the kcat in all the reactions of the model.

% Get the kcat value based on all the kcats in the model
standardKcat = median(model.ec.kcat, 'omitnan');

enzSubSystems = cell(numel(model.ec.rxns), 1);

% Get the subSystem for the rxns with a GPR
for i = 1:numel(model.ec.rxns)
    idx = strcmpi(model.rxns, model.ec.rxns{i});
    enzSubSystems(i,1) = model.subSystems{idx};
end

% Remove from the list those with kcat zero
rxnsKcatZero = model.ec.kcat > 0;

% Determine the subSystems in model.ec
[enzSubSystem_group, enzSubSystem_names] = findgroups(enzSubSystems(rxnsKcatZero));

% Calculate the mean kcat value for each subSystem in model.ec
kcatSubSystem = splitapply(@mean, model.ec.kcat(rxnsKcatZero), enzSubSystem_group);

% Calculate the number of reactions for each subSystem in model.ec
numRxnsSubSystem = splitapply(@numel, model.ec.rxns(rxnsKcatZero), enzSubSystem_group);

% Find subSystems which contains < threshold rxns and assign a standard value
lowerThanTresh = numRxnsSubSystem < threshold;
kcatSubSystem(lowerThanTresh) = standardKcat;

% Find reactions without GPR
rxnsMissingGPR = find(cellfun(@isempty, model.grRules));

% Get and remove exchange, transport, spontaneous and pseudo reactions
[~, exchangeRxns]  = getExchangeRxns(model);
transportRxns = getTransportRxns(model);
[spontaneousRxns, ~] = modelAdapter.getSpontaneousReactions(model);
pseudoRxns = contains(model.rxnNames,'pseudoreaction');

rxnsMissingGPR(ismember(rxnsMissingGPR, exchangeRxns)) = [];
rxnsMissingGPR(ismember(rxnsMissingGPR, find(transportRxns))) = [];
rxnsMissingGPR(ismember(rxnsMissingGPR, find(spontaneousRxns))) = [];
rxnsMissingGPR(ismember(rxnsMissingGPR, find(pseudoRxns))) = [];

% Add a new metabolite named prot_standard
proteinStdMets.mets         = 'prot_standard';
proteinStdMets.metNames     = proteinStdMets.mets;
proteinStdMets.compartments = 'c';
if isfield(model,'metNotes')
    proteinStdMets.metNotes     = 'Standard enzyme-usage pseudometabolite';
end

model = addMets(model, proteinStdMets);

% Add a new gene to be consistent with ec field named standard
proteinStdGenes.genes = 'standard';
if isfield(model,'geneShortNames')
    proteinStdGenes.geneShortNames = 'std';
end

model = addGenesRaven(model, proteinStdGenes);

% Update .ec structure in model
model.ec.genes(end+1)      = {'standard'};
model.ec.enzymes(end+1)    = {'standard'};
model.ec.mw(end+1)         = standardMW;
model.ec.sequence(end+1)   = {''};
% Additional info
if isfield(model.ec,'concs')
    model.ec.concs(end+1)  = nan();
end

% Expand the enzyme rxns matrix
model.ec.rxnEnzMat =  [model.ec.rxnEnzMat, zeros(length(model.ec.rxns), 1)]; % 1 new enzyme
model.ec.rxnEnzMat =  [model.ec.rxnEnzMat; zeros(length(rxnsMissingGPR), length(model.ec.enzymes))]; % new rxns

numRxns = length(model.ec.rxns);
stdMetIdx = find(strcmpi(model.ec.enzymes, 'standard'));

for i = 1:numel(rxnsMissingGPR)
    rxnIdx = rxnsMissingGPR(i);
    kcatSubSystemIdx =  strcmpi(enzSubSystem_names, model.subSystems{rxnIdx});

    % Update .ec structure in model
    model.ec.rxns(end+1)     = model.rxns(rxnIdx);
    if all(kcatSubSystemIdx)
        model.ec.kcat(end+1) = kcatSubSystem(kcatSubSystemIdx);
    else
        model.ec.kcat(end+1) = standardKcat;
    end
    model.ec.source(end+1)   = {'Standard kcat'};
    model.ec.notes(end+1)    = {''};
    model.ec.eccodes(end+1)  = {''};

    % Update the enzyme rxns matrix
    model.ec.rxnEnzMat(numRxns+i, stdMetIdx) = 1;
end
end
