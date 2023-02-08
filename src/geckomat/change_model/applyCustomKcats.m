function [model, rxnUpdated, rxnNotMatch] = applyCustomKcats(model, modelAdapter)
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
%   rxnUpdated        list of updated reactions, new kcats were applied 
%   rxnNotMatch      list of reaction which the custom information provided does not have 
%                             full match (100%) based on GPR rules, suggested to be curated by the
%                             user
%
% Usage:
%   [model, rxnUpdated, rxnNotMatch] = applyComplexData(model, modelAdapter);

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

                % Only keep reactions that have full match with the gene
                % rule
                [C,ia,ib] = intersect(enzIdx, allEnzInRxn);
                if numel(ia) == numel(enzIdx) && numel(ib) == numel(allEnzInRxn)
                    rxnToUpdate(rxnIdxs(j)) = 1;
                    model.ec.kcat(rxnIdxs(j)) = customKcats{4}(i);

                    % Add note mentioning manual kcat change
                    if  isempty(model.ec.notes{rxnIdxs(j), 1})
                        model.ec.notes{rxnIdxs(j), 1} = 'kcat modified manually';
                    else
                        model.ec.notes{rxnIdxs(j), 1} = [model.ec.notes{rxnIdxs(j), 1} ', kcat modified manually'];
                    end
                else
                    rxnNotMatch(rxnIdxs(j)) = 1;
                    disp(['Reaction ' model.ec.rxns{rxnIdxs(j)} ' does not have full match with the custom data: ' customKcats{1}{i}]);
                end
            end
        end
    end
    
    % Apply the new kcat values to the model
    model = applyKcatConstraints(model, rxnToUpdate);
    rxnUpdated = model.ec.rxns(find(rxnToUpdate));
    rxnNotMatch = model.ec.rxns(find(rxnNotMatch));

end

end

