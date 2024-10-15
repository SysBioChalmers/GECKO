function [model, rxnUpdated, notMatch] = applyCustomKcats(model, customKcats, modelAdapter)
% applyCustomKcats
%   Apply user defined kcats.  Reads data/customKcats.tsv in the obj.params.path
%   specified in the model adapter. Alternatively, a customKcats structure can
%   provided, as specified below.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   customKcats     structure with custom kcat information. If nothing
%                   is provided, an attempt will be made to read
%                   data/customKcats.tsv from the obj.params.path folder
%                   specified in the modelAdapter.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
% Output:
%   model           ecModel where kcats for defined proteins have been
%                   changed
%   rxnUpdated      ids list of updated reactions, new kcats were applied
%   notMatch        table with the list of reactions which the custom
%                   information provided does not have full match (> 50%)
%                   based on GPR rules. Then, they are suggested to be
%                   curated by the user
%
%   customKcats structure:
%   - proteins    protein identifiers, multiple for the same kcat (in case
%                 of a protein complex) are separated by ' + '
%   - genes       gene identifiers (optional, not used in matching)
%   - gene_name   short gene name (optional, not used in matching)
%   - kcat        new kcat value (one per entry)
%   - rxns        reaction identifiers, multiple for the same kcat are
%                 separated by ',' (see further explanation below)
%   - notes       will be appended to model.ec.notes (optional)
%   - stoicho     complex stoichiometry, separated by ' + ' (examples: '1'
%                 or '3 + 1'), matching the order in proteins field
%   
%   Matching order:
%   (1) reactions are identified by .proteins and .rxns
%   (2) reactions are identified by .proteins only (empty .rxns entry)
%       no additional checks are made: is a reaction annotated with these
%       proteins? => its kcat will be updated, irrespective of the exact
%       reaction, direction, substrate, etc.
%   (3) reactions are identified by .rxns only (empty .proteins entry)
%       no additional checks are made: is a reaction derived from the
%       original reaction identifier => its kcat will be updated,
%       irrespective of the annotated protein
%
%   customKcats.rxns field:
%   The reaction identifiers are from the ORIGINAL model, before _EXP_
%   suffixes were added by makeEcModel. Reaction directionality IS however
%   specified, with a _REV suffix.
%   Example entries:
%   'r_0001'     will match r_0001, r_0001_EXP_1, r_0001_EXP_2 etc., but
%                not r_0001_REV, r_0001_EXP_1_REV etc.
%   'r_0001_REV' will match r_0001_REV, r_0001_EXP_1_REV, r_0001_EXP_2_REV,
%                etc., but not r_0001, r_0001_EXP_1 etc.
%   Multiple identifiers should be comma separated (e.g. r_0001, r_0002)
%
% Usage:
%   [model, rxnUpdated, notMatch] = applyCustomKcats(model, customKcats, modelAdapter);

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

params = modelAdapter.params;

if nargin<2 || isempty(customKcats)
    customKcats = fullfile(params.path,'data','customKcats.tsv');
end
if isstruct(customKcats)
    if ~all(isfield(customKcats,{'proteins','kcat','rxns'}))
        error('The customKcats structure does not have all essential fields.');
    end
elseif isfile(customKcats)
    fID = fopen(customKcats, 'r');
    fileContent = textscan(fID, '%s %s %s %f %q %s %s', 'Delimiter', '\t', 'HeaderLines', 1);
    fclose(fID);
    clear customKcats
    customKcats.proteins    = fileContent{1};
    customKcats.genes       = fileContent{2};
    customKcats.gene_name   = fileContent{3};
    customKcats.kcat        = fileContent{4};
    customKcats.rxns        = fileContent{5};
    customKcats.notes       = fileContent{6};
    customKcats.stoicho     = fileContent{7};
else
    error(['Cannot find file: ' customKcats]);
end

rxnToUpdate = false(length(model.ec.rxns),1);
rxnNotMatch = false(length(model.ec.rxns),1);
evaluatedRule = cell(length(model.ec.rxns),1);
enzInModel = cell(length(model.ec.rxns),1);
ecRxnNoSuffix = regexprep(model.ec.rxns,'_EXP_\d+$','');

% Implementation for full GECKO formulation
if ~model.ec.geckoLight
    for i = 1:numel(customKcats.proteins)
        if isempty(customKcats.proteins{i})
            %If only reaction ID(s) is/are specified (and no proteins),
            %then apply the kcat to all isozymic reactions
            rxns    = strtrim(strsplit(customKcats.rxns{i}, ','));
            rxnIdxs = ismember(ecRxnNoSuffix,rxns);
            rxnToUpdate(rxnIdxs) = 1;
            model.ec.kcat(rxnIdxs) = customKcats.kcat(i);
        else
            prots = strtrim(strsplit(customKcats.proteins{i}, '+'));

            % Find the index for the enzymes which kcat will be changed. In
            % case the protein is not in the model generate a warning.
            try
                enzIdx = cellfun(@(x) find(strcmpi(model.ec.enzymes, x)), prots);
            catch
                enzIdx = [];
                printOrange(['WARNING: Protein(s) ' customKcats.proteins{i} ' were not found in the model.']);
            end

            % if not specific reactions are defined, find all the reaction
            % index where the enzyme is used
            if isempty(customKcats.rxns{i})
                temp_rxnIdxs = arrayfun(@(x) find(model.ec.rxnEnzMat(:, x)), enzIdx, 'UniformOutput', false);
                % otherwhise, If a set of reactions if defined, only get the index for those
                % but ignore any _EXP_ suffixes
            else
                rxns = strtrim(strsplit(customKcats.rxns{i}, ','));
                temp_rxnIdxs = arrayfun(@(x) find(strcmpi(ecRxnNoSuffix, x)), rxns, 'UniformOutput', false);
            end

            if ~isempty(temp_rxnIdxs)
                rxnIdxs = [];
                for j = 1:numel(temp_rxnIdxs)
                    rxnIdxs = [rxnIdxs; temp_rxnIdxs{j}];
                end

                % Check when multiple proteins are involved, since it can return same rxn n times
                rxnIdxs = unique(rxnIdxs); %unique(rxnIdxs{1, :});

                for j = 1:numel(rxnIdxs)
                    % Get all the enzymes involved in the reaction
                    allEnzInRxn = find(model.ec.rxnEnzMat(rxnIdxs(j),:));

                    C = intersect(enzIdx, allEnzInRxn);

                    % Determine the match percentage bewteen the rules
                    if numel(C) == numel(enzIdx) && numel(C) == numel(allEnzInRxn)
                        match = 1;
                    else
                        if numel(enzIdx) < numel(allEnzInRxn)
                            match = numel(C) / numel(allEnzInRxn);
                        else
                            match = numel(C) / numel(enzIdx);
                        end
                    end

                    % Update the kcat value stored in the model, if match 100%,
                    % otherwhise if >= 0.5 inform to the user to be validated
                    if match == 1
                        rxnToUpdate(rxnIdxs(j)) = 1;
                        model.ec.kcat(rxnIdxs(j)) = customKcats.kcat(i);

                        % Add note mentioning manual kcat change
                        model.ec.source{rxnIdxs(j),1} = 'custom';
                        if isfield(customKcats,'notes')
                            if isempty(model.ec.notes{rxnIdxs(j), 1}) && ~isempty(customKcats.notes{i})
                                model.ec.notes{rxnIdxs(j), 1} = customKcats.notes{i};
                            else
                                model.ec.notes{rxnIdxs(j), 1} = [model.ec.notes{rxnIdxs(j), 1} ', ' customKcats.notes{i}];
                            end
                        end
                    elseif match >= 0.5 && match < 1
                        rxnNotMatch(rxnIdxs(j)) = 1;
                        evaluatedRule{rxnIdxs(j), 1} = customKcats.proteins{i};
                        enzInModel{rxnIdxs(j), 1} = strjoin(model.ec.enzymes(allEnzInRxn), ' + ');
                    end
                end
            end
        end
    end

    % Apply the new kcat values to the model
    if ~isempty(find(rxnToUpdate, 1))
        model = applyKcatConstraints(model, rxnToUpdate);
    else
        printOrange('WARNING: No matches found. Consider checking the IDs or proteins in customKcats.');
    end

    rxnUpdated = model.ec.rxns(find(rxnToUpdate));

    % Remove from notMatch those reactions already updated
    remove = and(rxnToUpdate, rxnNotMatch);
    rxnNotMatch(remove) = 0;
    evaluatedRule(remove) = '';
    enzInModel(remove) = '';

    idRxns = model.ec.rxns(find(rxnNotMatch));
    fullIdx = cellfun(@(x) find(strcmpi(model.rxns, x)), idRxns);
    rxnsNames = model.rxnNames(fullIdx);
    evaluatedRule = evaluatedRule(~cellfun('isempty', evaluatedRule));
    enzInModel = enzInModel(~cellfun('isempty', enzInModel));
    rules = model.grRules(fullIdx);
    notMatch = table(idRxns, rxnsNames, evaluatedRule, enzInModel, rules, ...
        'VariableNames',{'rxns', 'name', 'custom enzymes', 'enzymes in model', 'rules'});
end
end
