function [model, rxnsMissingGPR, standardMW, standardKcat, rxnsNoKcat] = getStandardKcat(model, modelAdapter, threshold, fillZeroKcat)
% getStandardKcat
%   Calculate an standard kcat and standard molecular weight (MW) that can
%   be used to apply enzyme constraints to reactions without any associated
%   enzymes. Such reactions have either an empty model.grRules field, or
%   they have no match in model.ec.rxns (which can be the case if the genes
%   in the model.grRules field could not be mapped to enzymes). This is
%   done by adding those reactions to model.ec, assign a "standard"
%   pseudoenzyme with the standard MW (median of all proteins in the
%   organism) and standard kcat (median from all kcat, or subsystem
%   specific kcat).
%
%   A reaction is assigned a subSystem specific kcat values if the model
%   has a subSystems field and the reaction is annotated with a subSystem.
%   Only the first subSystem will be considered if multiple are annotated
%   to the same reaction.
%
%   Exchange, transport and pseudoreactions are filtered out, plus any
%   reaction identifiers specified in /data/pseudoRxns.tsv in the model
%   adapter folder.
%
%   In addition, reactions that are annotated with an enzyme (and therefore
%   already in model.ec), but not assigned any reaction-specific kcat value
%   (their model.ec.kcat entry is either 0 or NaN), can be assigned
%   standard kcat values by a similar approach. However, those reactions
%   will not be linked to the "standard" pseudoenzyme, but will use the
%   enzyme that they had already been associated with.
%
%   Any pre-existing standard kcat assignments (identified by 'standard'
%   entires in model.ec.source) are removed when applying this function.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use
%                   the default model adapter).
%   threshold       a threshold to determine when use a kcat value based on
%                   the mean kcat of the reactions in the same subSystem or
%                   based on the median value of all the kcat in the model.
%                   Second option is used when the number of reactions in a
%                   determined subSystem is < threshold. (Optional, default
%                   = 10)
%   fillZeroKcat    logical whether zero kcat values should be replaced
%                   with standard kcat values. (Optional, default = true).
%
% Output:
%   model           ecModel where model.ec is expanded with a standard
%                   protein with standard kcat and standard MW, assigned to
%                   reactions without gene associations.
%   rxnsMissingGPR  a list of updated rxns identifiers with a standard value
%   standardMW      the standard MW value calculated
%   standardKcat    the standard Kcat value calculated 
%   rxnsNoKcat      a list of rxns identifiers whose zero kcat has been replaced
%
%   While model.ec.kcat is populated, applyKcatConstraints would still need
%   to be run to apply the new constraints to the S-matrix.
%
% Usage:
%    [model, rxnsMissingGPR, standardMW, standardKcat, rxnsNoKcat] = getStandardKcat(model, modelAdapter, threshold, fillZeroKcat);

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 3 || isempty(threshold)
    threshold = 10;
end

if nargin < 4 || isempty(fillZeroKcat)
    fillZeroKcat = true;
end

databases = loadDatabases('uniprot', modelAdapter);

% An standard MW is defined for all the rxns which does not have a GPR
% rule defined. This is based in all the proteins reported for the specific
% organism in uniprot
standardMW = median(databases.uniprot.MW, 'omitnan');

% An standard Kcat is defined for all the rxns which does not have a GPR
% rule defined. In this case, the kcat value for a particular reaction is
% defined as the mean of the kcat values of the reactions involved in the
% same subsystem in which the given reaction is involved. Nevertheless, if
% a subSystem have a number of reactions lower than a treshold, the kcat
% value will be the median of the kcat in all the reactions of the model.

% Remove from the list those with kcat zero
rxnsKcatZero = model.ec.kcat > 0;

% Get the kcat value based on all the kcats in the model
standardKcat = median(model.ec.kcat(rxnsKcatZero), 'omitnan');

% If the model have subSystems assigned calculate kcat based on subSystem
if isfield(model,'subSystems') && ~all(cellfun(@isempty, model.subSystems))
    standard = false;
    if model.ec.geckoLight
        modRxns = extractAfter(model.ec.rxns,4);
    else
        modRxns = model.ec.rxns;
    end
    % Map ec-rxns to model.rxns
    [~,rxnIdx]  = ismember(modRxns,model.rxns);
    % Choose first subSystem
    enzSubSystems = flattenCell(model.subSystems(rxnIdx));
    enzSubSystems = enzSubSystems(:,1);
    if ~all(cellfun(@isempty, enzSubSystems))

    % Make list of unique subsystems, and which rxns are linked to them
    [enzSubSystem_names, ~, rxnToSub] = unique(enzSubSystems);
    % Make matrix of ec-rxns vs. unique subsystem index
    ind = sub2ind([numel(enzSubSystem_names) numel(enzSubSystems)],rxnToSub',1:numel(rxnToSub));
    kcatSubSystem = false([numel(enzSubSystem_names) numel(enzSubSystems)]);
    kcatSubSystem(ind) = true;
    % Number of kcats per subSystem
    kcatsPerSubSystem = sum(kcatSubSystem,2);
    % Calculate average kcat values per subSystem
    kcatSubSystem = (kcatSubSystem*model.ec.kcat)./kcatsPerSubSystem;
    kcatSubSystem(kcatsPerSubSystem < threshold) = standardKcat;
    else
        standard = true;
        printOrange('WARNING: No subSystem-specific kcat values can be calculated');
    end
else
    standard = true;
    printOrange('WARNING: No subSystem-specific kcat values can be calculated');
end

% Find reactions without GPR
rxnsMissingGPR = find(cellfun(@isempty, model.grRules));

% Find reactions with GPR but without model.ec entry (for instance due to
% no protein matching)
rxnsMissingEnzyme = find(~cellfun(@isempty, model.grRules));
if model.ec.geckoLight
    ecRxnsList = unique(extractAfter(model.ec.rxns,4));
else
    ecRxnsList = model.ec.rxns;
end
rxnsMissingEnzyme = find(and(~ismember(model.rxns(rxnsMissingEnzyme),ecRxnsList), ~contains(model.rxns(rxnsMissingEnzyme),'usage_prot_')));
rxnsMissingGPR = [rxnsMissingGPR;rxnsMissingEnzyme];

% Get custom list of reaction IDs to ignore, if existing. First column
% contains reaction IDs, second column contains reaction names for
% reference only.
if exist(fullfile(params.path,'data','pseudoRxns.tsv'),'file')
    fID        = fopen(fullfile(params.path,'data','pseudoRxns.tsv'));
    fileData   = textscan(fID,'%s %s','delimiter','\t');
    fclose(fID);
    customRxns = fileData{1};
    customRxns = find(ismember(model.rxns,customRxns));
else
    customRxns = [];
end
% Get and remove exchange, transport, spontaneous and pseudo reactions
[~, exchangeRxns]  = getExchangeRxns(model);
transportRxns = getTransportRxns(model);
[spontaneousRxns, ~] = modelAdapter.getSpontaneousReactions(model);
pseudoRxns = contains(model.rxnNames,'pseudoreaction');
slimeRxns = contains(model.rxnNames,'SLIME rxn');
rxnsToIgnore = [customRxns; exchangeRxns; find(transportRxns); ...
                find(spontaneousRxns); find(pseudoRxns); find(slimeRxns)];
rxnsMissingGPR(ismember(rxnsMissingGPR, rxnsToIgnore)) = [];

% Only add if not geckoLight & getStandardKcat was not run earlier
if ~any(strcmp(model.mets,'prot_standard'))
    % Add a new gene to be consistent with ec field named standard
    proteinStdGenes.genes = 'standard';
    if isfield(model,'geneShortNames')
        proteinStdGenes.geneShortNames = 'std';
    end
    model = addGenesRaven(model, proteinStdGenes);

    if ~model.ec.geckoLight
        % Add a new metabolite named prot_standard
        proteinStdMets.mets         = 'prot_standard';
        proteinStdMets.metNames     = proteinStdMets.mets;
        % Validate compartment
        compartmentID = strcmp(model.compNames,params.enzyme_comp);
        if ~any(compartmentID)
            error(['Compartment ' params.enzyme_comp ' (specified in params.enzyme_comp) '...
                'cannot be found in model.compNames'])
        end
        proteinStdMets.compartments = model.comps(compartmentID);

        if isfield(model,'metNotes')
            proteinStdMets.metNotes = 'Standard enzyme-usage pseudometabolite';
        end
        model = addMets(model, proteinStdMets);

        % Add a protein usage reaction if not a light version
        proteinStdUsageRxn.rxns         = {'usage_prot_standard'};
        proteinStdUsageRxn.rxnNames     = proteinStdUsageRxn.rxns;
        proteinStdUsageRxn.mets         = {proteinStdMets.mets, 'prot_pool'};
        proteinStdUsageRxn.stoichCoeffs = [-1, 1];
        proteinStdUsageRxn.lb           = -1000;
        proteinStdUsageRxn.ub           = 0;
        proteinStdUsageRxn.grRules      = proteinStdGenes.genes;

        model = addRxns(model, proteinStdUsageRxn);
    end
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
end
stdMetIdx = find(strcmpi(model.ec.enzymes, 'standard'));

% Remove previous standard kcat assignment
oldStandardEnz = find(strcmp(model.ec.source,'standard'));
if ~isempty(oldStandardEnz)
    oldStandardProt = logical(model.ec.rxnEnzMat(oldStandardEnz,stdMetIdx));
    % If annotated with real enzyme => set kcat to zero
    model.ec.kcat(oldStandardEnz(~oldStandardProt))        = 0;
    model.ec.source(oldStandardEnz(~oldStandardProt))      = {''};
    % If annotated with standard protein => remove entry
    model.ec.rxns(oldStandardEnz(oldStandardProt))        = [];
    model.ec.kcat(oldStandardEnz(oldStandardProt))        = [];
    model.ec.source(oldStandardEnz(oldStandardProt))      = [];
    model.ec.notes(oldStandardEnz(oldStandardProt))       = [];
    model.ec.eccodes(oldStandardEnz(oldStandardProt))     = [];
    model.ec.rxnEnzMat(oldStandardEnz(oldStandardProt),:) = [];
end

numRxns = length(model.ec.rxns);
for i = 1:numel(rxnsMissingGPR)
    rxnIdx = rxnsMissingGPR(i);

    % Update .ec structure in model
    if ~model.ec.geckoLight
        model.ec.rxns(end+1)     = model.rxns(rxnIdx);
        % Add prefix in case is light version
    else
        model.ec.rxns{end+1}     = ['001_' model.rxns{rxnIdx}];
    end

    if ~standard
        kcatSubSystemIdx = strcmpi(enzSubSystem_names, model.subSystems{rxnIdx}(1));
        if all(kcatSubSystemIdx)
            model.ec.kcat(end+1) = kcatSubSystem(kcatSubSystemIdx);
        else
            model.ec.kcat(end+1) = standardKcat;
        end
    else
        model.ec.kcat(end+1) = standardKcat;
    end

    model.ec.source(end+1)   = {'standard'};
    model.ec.notes(end+1)    = {''};
    model.ec.eccodes(end+1)  = {''};

    % Update the enzyme rxns matrix
    model.ec.rxnEnzMat(numRxns+i, stdMetIdx) = 1;
end
% Get the rxns identifiers of the updated rxns
rxnsMissingGPR = model.rxns(rxnsMissingGPR);

if fillZeroKcat
    zeroKcat = model.ec.kcat == 0 | isnan(model.ec.kcat);
    model.ec.kcat(zeroKcat)     = standardKcat;
    model.ec.source(zeroKcat)   = {'standard'};
    rxnsNoKcat = model.ec.rxns(zeroKcat);
else
    rxnsNoKcat = [];
end
end

function Cflat = flattenCell(C,strFlag)
%FLATTENCELL  Flatten a nested column cell array into a matrix cell array.
%
% CFLAT = FLATTENCELL(C) takes a column cell array in which one or more
% entries is a nested cell array, and flattens it into a 2D matrix cell
% array, where the nested entries are spread across new columns.
%
% CFLAT = FLATTENCELL(C,STRFLAG) if STRFLAG is TRUE, empty entries in the
% resulting CFLAT will be replaced with empty strings {''}. Default = FALSE
if nargin < 2
    strFlag = false;
end

% determine which entries are cells
cells = cellfun(@iscell,C);

% determine number of elements in each nested cell
cellsizes = cellfun(@numel,C);
cellsizes(~cells) = 1;  % ignore non-cell entries

% flatten single-entry cells
Cflat = C;
Cflat(cells & (cellsizes == 1)) = cellfun(@(x) x{1},Cflat(cells & (cellsizes == 1)),'UniformOutput',false);

% iterate through multi-entry cells
multiCells = find(cellsizes > 1);
for i = 1:length(multiCells)
    cellContents = Cflat{multiCells(i)};
    Cflat(multiCells(i),1:length(cellContents)) = cellContents;
end

% change empty elements to strings, if specified
if ( strFlag )
    Cflat(cellfun(@isempty,Cflat)) = {''};
end
end
