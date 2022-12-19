%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eccodes = getECCodes(rxnEnzMat, genes, modelAdapter, action)
% Retrieves the enzyme codes for each of the reactions for a given genome
% scale model (GEM). This function was adapted from the previous GetEnzymeCodes in GECKO 2.
%
% INPUT:    enzGeneMap      From the ec structure in the model
%           genes           From the ec structure in the model
%           modelAdapter    Model adapter (for example an instance of YeastGEMAdapter)
%           action:         Response action if multiple proteins with
%                           different EC numbers are found for a given gene in
%                           a metabolic reaction.
%                           - 'display' Displays all found multiplicities.
%                           - 'ignore'  Ignore multiplicities and use the
%                              protein with the lowest index in the database.
%                           - 'add'     Adds all the multiple proteins as
%                                       isoenzymes for the given reaction.
%           
% OUTPUTS:  eccodes - can be a single EC code or several separated by ';'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eccodes = getECCodes(rxnEnzMap, genes, modelAdapter, action)

if nargin<4
    action = 'display';
end

fprintf('Retrieving EC numbers...')

data      = load(modelAdapter.getFilePath('ProtDatabase.mat'));
swissprot = data.swissprot;
kegg      = data.kegg;

swissprot = standardizeDatabase(swissprot);
kegg      = standardizeDatabase(kegg);

DBprotSwissprot     = swissprot(:,1);
DBgenesSwissprot    = flattenCell(swissprot(:,3));
DBecNumSwissprot    = swissprot(:,4);
DBMWSwissprot       = swissprot(:,5);

DBprotKEGG          = kegg(:,1);
DBgenesKEGG         = flattenCell(kegg(:,3));
DBecNumKEGG         = kegg(:,4);
DBMWKEGG            = kegg(:,5);

n = size(rxnEnzMap,1);

eccodes = cell(n,1);

conflicts  = cell(1,4);

%Build an index from gene to prot for faster processing later
[x,y] = size(DBgenesSwissprot);
genesForIndex = reshape(DBgenesSwissprot, x*y, 1);
genesForIndex = genesForIndex(~cellfun(@isempty, genesForIndex));
genesForIndex = unique(genesForIndex); %18360
geneIndex = cell(length(genesForIndex),1);
geneHashMap = containers.Map(genesForIndex,1:length(genesForIndex));
protIndices = 1:length(DBgenesSwissprot(:,1));
for i = 1:y
    tmp1 = DBgenesSwissprot(:,i);
    sel = ~cellfun(@isempty, tmp1);
    indices = cell2mat(values(geneHashMap,tmp1(sel)));
    protIndicesSel = protIndices(sel);
    for j = 1:length(indices)
        geneIndex{indices(j)} = [geneIndex{indices(j)};protIndicesSel(j)]; 
    end
end

rxnEnzMap = logical(rxnEnzMap);

for i = 1:n
    gns = genes(rxnEnzMap(i,:).');
    if ~isempty(gns)
        %Find match in Swissprot:
        [new_EC,multGenes] = findECInDB(gns,DBprotSwissprot,DBgenesSwissprot,DBecNumSwissprot,DBMWSwissprot,geneIndex,geneHashMap);
        if ~isempty(new_EC)
            DBase    = 'swissprot';
            if ~isempty(multGenes{1})
                multGenes{3} = DBase;
            end
        else
            %Find match in KEGG (skipped optimizing this step) - may be good to fix this later to get rid of the findInDBOld function:
            [new_EC,multGenes] = findECInDB(gns,DBprotKEGG,DBgenesKEGG,DBecNumKEGG,DBMWKEGG,geneIndex,geneHashMap);
            if ~isempty(new_EC)
                DBase    = 'kegg';
                if ~isempty(multGenes{1})
                    multGenes{3} = DBase;
                end
            end
        end
        eccodes{i} = new_EC; 
        
        if ~isempty(multGenes{1})
            %Rxn index
            conflicts{1} = [conflicts{1};i];
            %Gene IDs
            conflicts{2} = [conflicts{2};multGenes{1}];
            %Indexes in DB
            conflicts{3} = [conflicts{3};multGenes{2}];
            %DB name
            conflicts{4} = [conflicts{4};multGenes{3}];
            
            %{ I don't understand the purpose of this, let's skip it for now
            %if strcmpi(action,'add')
            %    if strcmpi(DBase,'swissprot')
            %        [uni,EC,MW,Genes] = addMultipleMatches(uni,EC,MW,Genes,multGenes,swissprot);
            %    elseif strcmpi(DBase,'KEGG')
            %        [uni,EC,MW,Genes] = addMultipleMatches(uni,EC,MW,Genes,multGenes,kegg);
            %    end
            %end
            %}
        end
    end
    if rem(i,500) == 0 || i == n
        fprintf('.')
    end
end

%Display error message with the multiple gene-protein matches found
if strcmpi(action,'display') && ~isempty(conflicts{1})
    displayErrorMessage(conflicts,swissprot,kegg)
end

fprintf(' Done!\n')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function database = standardizeDatabase(database)
for i = 1:length(database)
    database{i,3} = strsplit(database{i,3},' ');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = union_string(cell_array)
%Receives any 1xn cell array and returns the union of all non empty
%elements as a string
nonempty = ~cellfun(@isempty,cell_array);
str      = strjoin(cell_array(nonempty)',' ');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I don't understand this part, skipping it for now
%{
function [uni,EC,MW,genes] = addMultipleMatches(uni,EC,MW,genes,conflicts,DB)
for i=1:length(conflicts{1})
    indexes = conflicts{2}{i};
    for j=2:length(indexes)
        indx  = indexes(j);
        uni   = [uni; DB{indx,1}];
        ECset = getECstring('',DB{indx,4});
        EC    = [EC; {ECset}];
        MW    = [MW; DB{indx,5}];
        genes = [genes; conflicts{1}{i}];
    end
end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayErrorMessage(conflicts,swissprot,kegg)
STR = '\n Some  genes with multiple associated proteins were found, please';
STR = [STR, ' revise case by case in the protDatabase.mat file:\n\n'];   
for i=1:length(conflicts{1})
    if strcmpi(conflicts{4}{i},swissprot)
        DB = swissprot;
    else
        DB = kegg;
    end
    proteins = DB(conflicts{3}{i},1);
    STR = [STR, '- gene: ' conflicts{2}{i} '  Proteins: ' strjoin(proteins) '\n'];
end
STR = [STR, '\nIf a wrongly annotated case was found then call the '];
STR = [STR, 'getEnzymeCodes.m function again with the option action'];
STR = [STR, '= ignore\n\n'];
STR = [STR, 'If the conflicting proteins are desired to be kept as isoenzymes'];
STR = [STR, ' then call the getEnzymeCodes.m function'];
STR = [STR, ' again with the option action = add\n'];
error(sprintf(STR))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%