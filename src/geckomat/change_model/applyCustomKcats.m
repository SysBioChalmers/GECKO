function [model, rxnUpdated, notMatch] = applyCustomKcats(model, modelAdapter)
% applyCustomKcats
%   Apply user defined kcats
%
% Input:
%   model                an ecModel in GECKO 3 version
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                             default model adapter).
%
% Output:
%   model                ecModel where kcats for defined proteins have been
%                             changed
%   rxnUpdated        ids list of updated reactions, new kcats were applied
%   notMatch            table with the list of reactions which the custom information provided
%                             does not have full match (> 50%) based on GPR rules. Then, they are 
%                             suggested to be curated by the user
%
% Usage:
%   [model, rxnUpdated, notMatch] = applyComplexData(model, modelAdapter);

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

fID = fopen(fullfile(params.path, 'data', 'customKcats.tsv'), 'r');
customKcats = textscan(fID, '%s %s %s %f %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1);
fclose(fID);

rxnToUpdate = false(length(model.ec.rxns),1);
rxnNotMatch = false(length(model.ec.rxns),1);
evaluatedRule = cell(length(model.ec.rxns),1);
enzInModel = cell(length(model.ec.rxns),1);

% Implementation for full GECKO formulation
if ~model.ec.geckoLight
    for i = 1:numel(customKcats{1})
        prots = strtrim(strsplit(customKcats{1}{i}, '+'));

        % Find the index for the enzymes which kcat will be changed. In
        % case the protein is not in the model generate a warning.
        try
            enzIdx = cellfun(@(x) find(strcmpi(model.ec.enzymes, x)), prots);
        catch
            enzIdx = [];
            disp( ['Protein(s) ' customKcats{1}{i} ' were not found in the model.']);
        end

        % Find all the reactions index where the enzyme is used
        temp_rxnIdxs = arrayfun(@(x) find(model.ec.rxnEnzMat(:, x)), enzIdx, 'UniformOutput', false);

        if ~isempty(temp_rxnIdxs)
            rxnIdxs = [];
            for j = 1:numel(temp_rxnIdxs)
                rxnIdxs = [rxnIdxs; temp_rxnIdxs{j}];
            end
            %rxnIdxs = cell2mat(temp_rxnIdxs);%arrayfun(@(x) horzcat(rxnIdxs, x), temp_rxnIdxs);
            % Check when multiple proteins are involved, since it can return same rxn n times
            rxnIdxs = unique(rxnIdxs); %unique(rxnIdxs{1, :});

            for j = 1:numel(rxnIdxs)
                % Get all the enzymes involved in the reaction
                allEnzInRxn = find(model.ec.rxnEnzMat(rxnIdxs(j),:));

                C = intersect(enzIdx, allEnzInRxn);

                % Determine the match percentaje bewteen the rules
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
                    model.ec.kcat(rxnIdxs(j)) = customKcats{4}(i);

                    % Add note mentioning manual kcat change
                    if  isempty(model.ec.notes{rxnIdxs(j), 1})
                        model.ec.notes{rxnIdxs(j), 1} = 'kcat modified manually';
                    else
                        model.ec.notes{rxnIdxs(j), 1} = [model.ec.notes{rxnIdxs(j), 1} ', kcat modified manually'];
                    end
                elseif match >= 0.5 && match < 1
                    rxnNotMatch(rxnIdxs(j)) = 1;
                    evaluatedRule{rxnIdxs(j), 1} = customKcats{1}{i};
                    enzInModel{rxnIdxs(j), 1} = strjoin(model.ec.enzymes(allEnzInRxn), ' + ');
                end
            end
        end
    end

    % Apply the new kcat values to the model
    model = applyKcatConstraints(model, rxnToUpdate);

    rxnUpdated = model.ec.rxns(find(rxnToUpdate));

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

