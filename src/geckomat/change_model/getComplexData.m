function complexInfo = getComplexData(taxonomicID, modelAdapter)
% getComplexData
%   Download curated complex stochiometries from the EMBL-EBI Complex
%   Portal database. Writes data/ComplexPortal.json in the obj.params.path
%   specified in the model adapter.
%
% Input:
%   taxonomicID     taxonomic identifier for which complex data should be
%                   downloaded. Only taxonomic identifiers allowed are
%                   those included on Complex Portal:
%                   https://www.ebi.ac.uk/complexportal/complex/organisms
%                   If empty, no complex data is downloaded, if 0 (zero),
%                   complex data from all organisms in Complex Portal is
%                   downloaded.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
% Output:
%   complexInfo     structure with data downloaded from Complex Portal.
%                   Contains the following fields:
%                   - complexID: id of the complex on Complex Portal
%                   - name: name of the complex on Complex Portal
%                   - species: organism containing the complex
%                   - geneName: names of the genes in the complex
%                   - protID: Uniprot IDs of the proteins in the complex
%                   - stochiometry: the complex stochiometry given in the same
%                     order as the genes and proteins in geneName and protID
%                   - defined:  0 if Complex Portal has no defined stochiometry
%                               1 if defined subunit stochiometry is given
%                               2 if complex consists of sub-complexes, whose
%                                 subunit stochiometries are given
% Usage
%   complexInfo = getComplexData(organism, modelAdapter);

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if nargin<1 || isempty(taxonomicID)
    taxonomicID = modelAdapter.getParameters().complex.taxonomicID;
end

params = modelAdapter.params;
if isempty(taxonomicID) % Can be empty when gathered from model adapter
    printOrange('WARNING: No taxonomicID specified.')
    return
elseif taxonomicID == 0
    taxonomicID = [];
end

webOptions = weboptions('Timeout', 30);
try
    url1 = 'https://www.ebi.ac.uk/intact/complex-ws/search/*';
    if ~isempty(taxonomicID)
        url1 = [url1 '?facets=species&filters=species:("' num2str(taxonomicID) '")'];
    end
    data = webread(url1,webOptions);
catch ME
    if (strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError'))
        error('Cannot connect to the Complex Portal, perhaps the server is not responding');
    end
end
if data.size == 0
    error('No data could be gathered from Complex Portal for the specified taxonomicID.')
end
complexData = cell(data.size,7);

progressbar('Retrieving information for complexes');
for i = 1:data.size
    progressbar(i/data.size) % Update progress bar
    url2 = 'https://www.ebi.ac.uk/intact/complex-ws/complex/';
    complexID = data.elements(i,1).complexAC;
    try
        temp = webread([url2 complexID],webOptions);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError'))
            printOrange(['WARNING: Cannot retrieve the information for ' complexID '.\n']);
        end
        temp = [];
    end

    if ~isempty(temp)
        complexData(i,1) = {temp.complexAc};
        complexData(i,2) = {temp.name};
        complexData(i,3) = {temp.species};

        idxIntType = find(strcmpi({temp.participants.interactorType}, 'protein'));

        % Some complex reported are 'stable complex', then, save the id
        % complex and but set the genes and protein to a empty string.
        if numel(idxIntType) > 0
            complexData(i,4) = {{temp.participants(idxIntType).name}};
            complexData(i,5) = {{temp.participants(idxIntType).identifier}};
        else
            complexData(i,4) = {{temp.participants.name}};
            complexData(i,5) = {{temp.participants.identifier}};
        end

        % Portal complex has two stochiometry values, a minimum and
        % maximum value. Only minimum will be store. In some cases,
        % some complex does not have stochiometry coefficient, then, it
        % will be fill with zeros
        if ~cellfun('isempty',{temp.participants.stochiometry})
            % For some reason, if there is only one protein in the complex
            % split function does nor return a cell nx2, instead is 2x1,
            % then assign an incorrect stochiometry
            switch numel(idxIntType)
                case 0 % Contains complexes
                    stochiometry = split({temp.participants.stochiometry}.', ',');
                    complexData(i,7) = {2};
                case 1 % Contains one protein
                    stochiometry = split({temp.participants(idxIntType).stochiometry}.', ',').';
                    complexData(i,7) = {1};
                otherwise
                    stochiometry = split({temp.participants(idxIntType).stochiometry}.', ',');
                    complexData(i,7) = {1};
            end
            values = str2double(erase(stochiometry(:,1),"minValue: ")).';
            complexData(i,6) = {values};
        else
            complexData(i,6) = {repelem(0,numel(complexData{i,4}))};
            complexData(i,7) = {0};
        end
    end
end
fprintf('\n');

% Expand complexes of complexes
complexComplex = find([complexData{:,7}]==2);
if ~isempty(complexComplex)
    for i=1:numel(complexComplex)
        subComplex    = complexData{complexComplex(i),5};
        subComplexS   = complexData{complexComplex(i),6};
        subComplexIdx = find(ismember(complexData(:,1),subComplex));
        allGenes = horzcat(complexData{subComplexIdx,4});
        allProts = horzcat(complexData{subComplexIdx,5});
        allStoch = {complexData{subComplexIdx,6}};
        for j=1:numel(subComplex)
            allStoch{j}=allStoch{j}*subComplexS(j);
        end
        allStoch = horzcat(allStoch{:});
        [allGenes,ia,ic] = unique(allGenes,'stable');
        allProts = allProts(ia);
        allStoch = splitapply(@sum, allStoch', ic);
        complexData{complexComplex(i),4} = allGenes;
        complexData{complexComplex(i),5} = allProts;
        complexData{complexComplex(i),6} = allStoch;
    end
end

rowHeadings = {'complexID','name','species','geneName','protID','stochiometry','defined'};

complexInfo = cell2struct(complexData, rowHeadings, 2);

% Convert to a JSON file
jsontxt = jsonencode(cell2table(complexData, 'VariableNames', rowHeadings));
% Write to a JSON file
fid = fopen(fullfile(params.path,'data','ComplexPortal.json'), 'w');
fprintf(fid, '%s', jsontxt);
fclose(fid);
fprintf('Model-specific ComplexPortal database stored at %s\n',fullfile(params.path,'data','ComplexPortal.json'));
end
