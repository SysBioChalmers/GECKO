function complexInfo = getComplexData(organism, writeFile)

% getComplexData
%
%   Download portal complex database with curated stochiomtery
%   
%   Input
%       - organism: organism available in complex portal. e.g.
%           * 'all'. For all the data in complex portal
%           * 'Homo sapiens'
%           * 'Mus musculus'
%           * 'Saccharomyces cerevisiae'
%         Note: If the organism is not defined, then download all
%               available information.
%
%       - writeFile: define if save the data as JSON format file.
%         default = false
%
%   Output
%       - complexInfo: data downloaded from complex portal
%
%   Usage
%         complexInfo = getComplexData('Saccharomyces cerevisiae',true)

if isequal(organism, 'all')
    organism = [];
end

if nargin < 2
    writeFile = false;
end

complexInfo = [];

% Switch to a valid name in complex portal 
if ~isempty(organism)
    if strcmp(organism,'Saccharomyces cerevisiae')
        organism = 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)';
    elseif strcmp(organism,'Schizosaccharomyces pombe')
        organism = 'Schizosaccharomyces pombe (strain 972 / ATCC 24843)';
    elseif strcmp(organism,'Escherichia coli')
        organism = 'Escherichia coli (strain K12)';
    end 
end

try
    url1 = 'https://www.ebi.ac.uk/intact/complex-ws/search/*';
    if ~isempty(organism)
        url1 = [url1 '?facets=species_f&filters=species_f:("' organism '")'];
    end

    data = webread(url1);

catch ME
    if (strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError'))
        disp('Cannot connect to the complex portal, perhaps the server is not responding');
    end
end

complexData = cell(data.size,6);
for i = 1:data.size

    url2 = 'https://www.ebi.ac.uk/intact/complex-ws/complex/';
    complexID = data.elements(i,1).complexAC;

    disp(['Retrieving information for ' complexID '. Complex ' int2str(i) ' of ' int2str(data.size)]);

    try
        temp = webread([url2 complexID]);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError'))
            warning(['Cannot retrieve the information for ' complexID]);
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
            complexData(i,4) = {' '};
            complexData(i,5) = {' '};
        end

        % Portal complex has two stochiometry values, a minimum and
        % maximum value. Only minimum will be store. In some cases,
        % some complex does not have stochiometry coefficient, then, it
        % will be fill with zeros
        if ~cellfun('isempty',{temp.participants.stochiometry})
            % For some reason, if there is only one protein in the complex
            % split function does nor return a cell nx2, instead is 2x1,
            % then assign an incorrect stochiometry 
            if numel(idxIntType) > 1
                stochiometry = split({temp.participants(idxIntType).stochiometry}.', ',');
            else
                stochiometry = split({temp.participants(idxIntType).stochiometry}.', ',').';
            end 
            values = str2double(erase(stochiometry(:,1),"minValue: ")).';
            complexData(i,6) = {values};
        else
            complexData(i,6) = {repelem(0,numel(complexData{i,4}))};
        end
    else
        %
        complexData(i,1) = {complexID};
        complexData(i,2) = {' '};
        complexData(i,3) = {' '};
        complexData(i,4) = {' '};
        complexData(i,5) = {' '};
        complexData(i,6) = {0};
    end
end

complexData = sortrows(complexData,3);

rowHeadings = {'complexID','name','specie','geneName','protID','stochiometry'};

complexInfo = cell2struct(complexData, rowHeadings, 2);

if writeFile == true
    % Convert to a JSON file
    jsontxt = jsonencode(cell2table(complexData, 'VariableNames', rowHeadings));

    % Write to a JSON file
    fid = fopen('../../databases/complex_data.json', 'w');
    fprintf(fid, '%s', jsontxt);
    fclose(fid);
end

end
