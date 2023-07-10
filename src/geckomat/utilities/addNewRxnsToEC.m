function [model, rxnsAdded] = addNewRxnsToEC(model, newRxns, newEnzymes, modelAdapter)
% getEngineeredModel
%   Add new reaction to an enzyme-constrained model
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   newRxns         structure with the new reaction information as follow:
%            rxns               cell array with unique strings that
%                               identifies each reaction
%            rxnNames           cell array with the names of each reaction
%            equations          cell array with equation strings. Decimal
%                               coefficients are expressed as "1.2".
%                               Reversibility is indicated by "<=>" or
%                               "=>".
%            grRules            cell array with the gene-reaction
%                               relationship for each reaction. For example
%                               "(A and B) or (C)" means that the reaction
%                               could be catalyzed by a complex between
%                               A & B or by C on its own. All the genes
%                               have to be present in model.genes. Add
%                               genes with addGenesRaven before calling
%                               this function if needed (opt, default '')
%   newEnzymes      structure with the new enzymes information as follow:
%            enzymes            cell array with uniprot id
%            genes              cell array with the respective gene
%            mw                 cell array with the MW
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
%   NOTE: Write equations as follows:
%               'xylitol[c] + NAD[c] => D-xylulose[c] + NADH[c] + H+[c]'
%         Note that the reversibility defined in the equation will be used
%         to split to construct irreversible reactions (add _REV) and the
%         rules defined will be used to expand the model (add _EXP_n) 
%
% Output:
%   model           ecModel whit new reactions
%   rxnsAdded       cell array with the reactions added
%
% Usage:
%   [model, rxnsAdded] = getEngineeredModel(ecModel, newRxns, newEnzymes, modelAdapter);

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

% Validate if GECKO version 3
if ~isfield(model, 'ec')
    error('Input model is not an ecModel in GECKO 3 format (with ecModel.ec structure)')
end

if ~model.ec.geckoLight

    % 1. Add new Enzyme-usage pseudometabolite

    % Validate compartment to add protein pseudoreactions
    compartmentID = strcmp(model.compNames,params.enzyme_comp);
    if ~any(compartmentID)
        error(['Compartment ' params.enzyme_comp ' (specified in params.enzyme_comp) '...
            'cannot be found in model.compNames'])
    end
    compartmentID = model.comps(compartmentID);

    proteinMets.mets         = strcat('prot_', newEnzymes.enzymes);
    proteinMets.metNames     = proteinMets.mets;
    proteinMets.compartments = compartmentID;
    if isfield(model,'metMiriams')
        proteinMets.metMiriams = repmat({struct('name',{{'sbo'}},'value',{{'SBO:0000252'}})},numel(proteinMets.mets),1);
    end
    if isfield(model,'metCharges')
        proteinMets.metCharges = zeros(numel(proteinMets.mets),1);
    end
    proteinMets.metNotes = repmat({'Enzyme-usage pseudometabolite'},numel(proteinMets.mets),1);
    model = addMets(model,proteinMets);

    % 2. Add new genes
    genesToAdd.genes = newEnzymes.genes;
    genesToAdd.geneShortNames = genesToAdd.genes;
    model = addGenesRaven(model, genesToAdd);

    % 3. Add new reactions

    % Get reversible reactions to split
    toIrrev = find(cellfun(@(x) contains(x,'<=>'),newRxns.equations));
    if ~isempty(toIrrev)
        for i=1:numel(toIrrev)
            idx = toIrrev(i);
            newEquation = split(newRxns.equations(idx),'<=>');
            % Remove directionality
            newRxns.equations(idx) = erase(newRxns.equations(idx), '<');

            % Add a new reaction to the list to be added
            newRxns.rxns{end+1} = [newRxns.rxns{idx} '_REV'];
            newRxns.rxnNames{end+1} = [newRxns.rxnNames{idx} ' (reversible)'];
            newRxns.equations{end+1} = [strtrim(newEquation{2, 1}) ' => ' strtrim(newEquation{1, 1})];
            newRxns.grRules{end+1} = newRxns.grRules{idx};
        end
    end

    % Get reactions to separate isozymes
    toSplit = find(cellfun(@(x) contains(x,'or'),newRxns.grRules));
    if ~isempty(toSplit)
        for i=1:numel(toSplit)
            idx = toSplit(i);
            newRules = split(newRxns.grRules(idx),'or');
            newRules = strtrim(erase(newRules, {'(' ')'}));

            for j=1:numel(newRules)
                newRxns.rxns{end+1} = [newRxns.rxns{idx} '_EXP_' num2str(j)];
                newRxns.rxnNames(end+1) = newRxns.rxnNames(idx);
                newRxns.equations(end+1) = newRxns.equations(idx);
                newRxns.grRules(end+1) = newRules(j);
            end

            % Remove the original one
            newRxns.rxns(idx) = [];
            newRxns.rxnNames(idx) = [];
            newRxns.equations(idx) = [];
            newRxns.grRules(idx) = [];
        end
    end

    newRxns.lb = zeros(length(newRxns.rxns),1);
    newRxns.ub = zeros(length(newRxns.rxns),1) + 1000;

    model = addRxns(model, newRxns, 3, '', true, true);

    % 4. Add protein usage reactions.
    usageRxns.rxns            = strcat('usage_',proteinMets.mets);
    usageRxns.rxnNames        = usageRxns.rxns;
    usageRxns.mets            = cell(numel(usageRxns.rxns),1);
    usageRxns.stoichCoeffs    = cell(numel(usageRxns.rxns),1);
    for i=1:numel(usageRxns.mets)
        usageRxns.mets{i}         = {proteinMets.mets{i}, 'prot_pool'};
        usageRxns.stoichCoeffs{i} = [-1,1];
    end
    usageRxns.lb              = zeros(numel(usageRxns.rxns),1) - 1000;
    usageRxns.ub              = zeros(numel(usageRxns.rxns),1);
    usageRxns.rev             = ones(numel(usageRxns.rxns),1);
    usageRxns.grRules         = newEnzymes.genes;
    model = addRxns(model,usageRxns);

    % 5. Update .ec field with enzyme data
    model.ec.enzymes    = [model.ec.enzymes; newEnzymes.enzymes(:)];
    model.ec.genes      = [model.ec.genes; newEnzymes.genes(:)];
    model.ec.mw         = [model.ec.mw; newEnzymes.mw(:)];
    model.ec.sequence   = [model.ec.sequence; repmat({''},numel(newEnzymes.enzymes),1)]; 
    if isfield(model.ec,'concs')
        model.ec.concs  = [model.ec.concs; NaN(numel(newEnzymes.enzymes),1)];
    end

    % 6. Expand the enzyme rxns matrix
    rxnWithRule = find(~cellfun('isempty',newRxns.grRules));

    model.ec.rxnEnzMat =  [model.ec.rxnEnzMat, zeros(length(model.ec.rxns), length(newEnzymes.enzymes))]; % add n new enzymes
    model.ec.rxnEnzMat =  [model.ec.rxnEnzMat; zeros(length(rxnWithRule), length(model.ec.enzymes))]; % add n new rxns

    % 7. Update .ec field with rxn data

    % Get reactions to add 
    numRxns = length(model.ec.rxns);
    for i = 1:numel(rxnWithRule)
        rxnIdx = rxnWithRule(i);

        model.ec.rxns(end+1)    = newRxns.rxns(rxnIdx);
        model.ec.kcat(end+1)    = 0;
        model.ec.source(end+1)  = {''};
        model.ec.notes(end+1)   = {''};
        model.ec.eccodes(end+1) = {''};

        % Get rules
        genes = newRxns.grRules(rxnIdx);
        genes = strtrim(erase(split(genes,'and'), {'(' ')'}));
        [~,enzIdx] = ismember(genes, model.ec.genes);

        % Update the enzyme rxns matrix
        model.ec.rxnEnzMat(numRxns+i, enzIdx) = 1;
    end
    
    rxnsAdded = newRxns.rxns(:);
end

end
