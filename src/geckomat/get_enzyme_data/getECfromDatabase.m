function model = getECfromDatabase(model, action, ecRxns, modelAdapter)
% getECfromDatabase
%   Populates the model.ec.eccodes field with enzyme codes that are
%   extracted from UniProt and KEGG databases, as assigned to the proteins
%   that catalyze the specific reactions.
%
% Input:
%   model           ec-model in GECKO 3 format
%   action          response action if multiple proteins with different EC
%                   numbers are found for a given gene in a metabolic
%                   reaction (optional, default 'display')
%                   - 'display' displays all found multiplicities
%                   - 'ignore'  ignore multiplicities and use the protein
%                               with the lowest index in the database.
%                   - 'add'     adds all the multiple proteins as
%                               isoenzymes for the given reaction
%   ecRxns          logical of length model.ec.rxns that specifies for
%                   which reactions the existing model.ec.eccodes entry
%                   should be kept and not modified by this function
%                   (optional, by default all model.ec.eccodes entries
%                   are populated by this function)
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
% Output:
%   model           ec-model with populated model.ec.eccodes

if nargin < 2
    action = 'display';
end

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

rxnEnzMat = model.ec.rxnEnzMat;
genes = modelAdapter.getUniprotCompatibleGenes(model.ec.genes);

data    = loadDatabases('both', modelAdapter);
uniprot = data.uniprot;
kegg    = data.kegg;

DBgenesUniprot  = data.uniprot.genes;
DBecNumUniprot  = data.uniprot.eccodes;
DBMWUniprot     = data.uniprot.MW;
%Build an index from gene to prot for faster processing later
[geneIndexUniprot,geneHashMapUniprot] = hashGeneToProt(DBgenesUniprot);

if ~isempty(kegg)
    DBgenesKEGG     = data.kegg.genes;
    DBecNumKEGG     = data.kegg.eccodes;
    DBMWKEGG        = data.kegg.MW;
    [geneIndexKEGG,geneHashMapKEGG]       = hashGeneToProt(DBgenesKEGG);
end
n = size(rxnEnzMat,1);

eccodes   = cell(n,1);
eccodes(:)= {''};
conflicts = cell(1,4);

rxnEnzMat = logical(rxnEnzMat);

for i = 1:n
    gns = genes(rxnEnzMat(i,:).');
    if ~isempty(gns)
        %Find match in Uniprot:
        [new_EC,multGenes] = findECInDB(gns,DBecNumUniprot,DBMWUniprot,geneIndexUniprot,geneHashMapUniprot);
        if ~isempty(new_EC)
            DBase    = 'uniprot';
            if ~isempty(multGenes{1})
                multGenes{3} = DBase;
            end
        end
        if ~isempty(kegg) && (isempty(new_EC) || endsWith(new_EC,'-'))
            %Find match in KEGG
            [new_EC_kegg,multGenes] = findECInDB(gns,DBecNumKEGG,DBMWKEGG,geneIndexKEGG,geneHashMapKEGG);
            if ~isempty(new_EC_kegg)
                DBase    = 'kegg';
                if ~isempty(multGenes{1})
                    multGenes{3} = DBase;
                end
                new_EC=new_EC_kegg;
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
            conflicts{4} = [conflicts{4};{multGenes{3}}];

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
end

%Display error message with the multiple gene-protein matches found
if strcmpi(action,'display') && ~isempty(conflicts{1})
    displayErrorMessage(conflicts,uniprot,kegg)
end

if nargin < 3 || isempty(ecRxns) || all(ecRxns)
    model.ec.eccodes = eccodes;
else
    if ~isfield(model.ec,'eccodes')
        model.ec.eccodes(1:numel(model.ec.rxns),1) = {''};
    end
    %Probably faster to subset with ecRxns in the beginning of the script,
    %but this was at the moment simpler to implement.
    model.ec.eccodes(ecRxns) = eccodes(ecRxns);
end

function displayErrorMessage(conflicts,uniprot,kegg)
STR = ['\n ' num2str(length(conflicts{1})) ' genes with multiple associated proteins were found, please'];
STR = [STR, ' revise case by case in the uniprot and kegg files:\n\n'];
for i=1:length(conflicts{1})
    if strcmpi(conflicts{4}{i},'uniprot')
        DB = uniprot.ID;
    else
        DB = kegg.uniprot;
    end
    proteins = DB(conflicts{3}{i});
    STR = [STR, '- gene: ' conflicts{2}{i} '  Proteins: ' strjoin(proteins) '\n'];
end
STR = [STR, '\nIf a wrongly annotated case was found then call the '];
STR = [STR, 'getECfromDatabase.m function again with the option action'];
STR = [STR, '= ignore\n\n'];
STR = [STR, 'If the conflicting proteins are desired to be kept as isoenzymes'];
STR = [STR, ' then call the getECfromDatabase.m function'];
STR = [STR, ' again with the option action = add\n'];
error(sprintf(STR))
end

function [geneIndex,geneHashMap]=hashGeneToProt(proteinDB)

[x,y] = size(proteinDB);
genesForIndex = reshape(proteinDB, x*y, 1);
genesForIndex = genesForIndex(~cellfun(@isempty, genesForIndex));
genesForIndex = unique(genesForIndex);
geneIndex = cell(length(genesForIndex),1);
geneHashMap = containers.Map(genesForIndex,1:length(genesForIndex));
protIndices = 1:length(proteinDB(:,1));
for i = 1:y
    tmp1 = proteinDB(:,i);
    sel = ~cellfun(@isempty, tmp1);
    indices = cell2mat(values(geneHashMap,tmp1(sel)));
    protIndicesSel = protIndices(sel);
    for j = 1:length(indices)
        geneIndex{indices(j)} = [geneIndex{indices(j)};protIndicesSel(j)];
    end
end
end
end
