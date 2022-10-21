function [model, foundComplex, proposedComplex] = applyComplexData(model, fileName)

% applyComplexData
%
%   Apply stochiometry for complex in an ecModel
%   
%   Input
%       - model: an ecModel in GECKO 3 version
%
%       - fileName: the directory to a JSON format file containing the
%         information of the complex for the specific microorganism. The
%         file must contain:
%           * complexID: id of the complex
%           * name: the name of the complex
%           * geneName: a cell array with the gene names in the complex
%           * protID: a cell array with the protein ids in the complex
%           * stochiometry: a cell array with the stochiometry for each
%             protein in the complex
%
%         Note: complex data can be downloaded using the getComplexData
%         function.
%
%   Output
%       - model: the ecModel with the rxnEnzMat populated
%
%       - foundComplex: the list of complex matched between the model and
%         the complex data.
%
%       - proposedComplex: the list of complex that match in more than 75%
%       between the model and the complex data or the gene information in
%       the model differs from the complex data
%
%   Usage
%         [model_test, foundComplex, proposedComplex] = applyComplexData(ecModel, 'complex_data.json');

jsonStr = fileread(fileName);
complexData = jsondecode(jsonStr);

foundComplex = cell(0,7);
proposedComplex = cell(0,8);

for i = 1:numel(model.ec.rxns)

    idxRxn = find(strcmpi(model.rxns, model.ec.rxns{i}));

    genes = split(model.grRules(idxRxn), 'and');
    genes = strtrim(genes);

    bestComplexIdx = 0;
    bestMatch = 0;
    protIdsComplex = [];
    protIdsModel = [];

    if numel(genes) > 1

        [~,~,idxProts] = intersect(genes,model.ec.genes, 'stable');

        if ~isempty(idxProts)
            protIdsModel = model.ec.enzymes(idxProts);

            for j = 1:height(complexData)
                
                protIdsComplex = complexData(j).protID;

                countMatch = ismember(protIdsModel,protIdsComplex);

                % Let´s use the complex data as reference
                match = sum(countMatch,'all')/numel(protIdsComplex);

                % Check if the protID match with the complex data based on match % higher than
                % 75%. Pick the highest match.
                if match >= 0.75 && match > bestMatch
                    bestComplexIdx = j;
                    bestMatch = match*100;
                    % In some cases all the model proteins are in the
                    % complex data, but they are less than those in the
                    % complex data.
                    if match == 1 && numel(protIdsModel) == numel(protIdsComplex)
                        break
                    end
                end
            end

            if bestMatch >= 75
                % Only get data with full match in the model and the
                % complex data. In some cases all the model proteins are in
                % the complex data, but they are less than those in complex data
                if bestMatch == 100 && numel(complexData(bestComplexIdx).protID) == numel(protIdsComplex)
                    foundComplex(end+1,1) = {model.ec.rxns{i}};
                    foundComplex(end,2) = {complexData(bestComplexIdx).complexID};
                    foundComplex(end,3) = {complexData(bestComplexIdx).name};
                    foundComplex(end,4) = {complexData(bestComplexIdx).geneName};
                    foundComplex(end,5) = {protIdsModel};
                    foundComplex(end,6) = {complexData(bestComplexIdx).protID};
                    foundComplex(end,7) = {complexData(bestComplexIdx).stochiometry};

                    % Some complex match in 100% but there is not stochiometry
                    % reported. In this case, assign a value of 1.
                    assignS = [complexData(bestComplexIdx).stochiometry];
                    assignS(assignS == 0) = 1;
                    model.ec.rxnEnzMat(i,idxProts) = assignS;
                else
                    proposedComplex(end+1,1) = {model.ec.rxns{i}};
                    proposedComplex(end,2) = {complexData(bestComplexIdx).complexID};
                    proposedComplex(end,3) = {complexData(bestComplexIdx).name};
                    proposedComplex(end,4) = {complexData(bestComplexIdx).geneName};
                    proposedComplex(end,5) = {protIdsModel};
                    proposedComplex(end,6) = {complexData(bestComplexIdx).protID};
                    proposedComplex(end,7) = {complexData(bestComplexIdx).stochiometry};
                    proposedComplex(end,8) = {bestMatch};
                end

            end
        end
    end
    
end

rowHeadings = {'rxn', 'complexID','name','genes','protID_model','protID_complex','stochiometry'};

foundComplex = cell2table(foundComplex, 'VariableNames', rowHeadings);

proposedComplex = cell2table(proposedComplex, 'VariableNames', [rowHeadings 'match']);

disp(['A total of ' int2str(numel(foundComplex(:,1))) ' complex have full match, and ' int2str(numel(proposedComplex(:,1))) ' proposed.'])

end