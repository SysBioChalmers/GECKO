%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_data = getEnzymeCodes(model,action)
% Retrieves the enzyme codes for each of the reactions for a given genome
% scale model (GEM).
%
% INPUT:    *model:      a GEM file (.mat format)
%           *action:     Response action if multiple proteins with
%                        different EC numbers are found for a given gene in
%                        a metabolic reaction.
%                           - 'display' Displays all found multiplicities.
%                           - 'ignore'  Ignore multiplicities and use the
%                              protein with the lowest index in the database.
%                           - 'add'     Adds all the multiple proteins as
%                                       isoenzymes for the given reaction.
%           
% OUTPUTS:  model_data, which contains:
%           *model:      The standardized GEM
%           *substrates: Substrates associated for each rxn
%           *products:   Products associated, when rxn is reversible
%           *uniprots:   All possible uniprot codes, for each rxn
%           *EC_numbers: All possible EC numbers, for each uniprot
%           *count(1):   #rxns with data from Swissprot
%           *count(2):   #rxns with data from KEGG
%           *count(3):   #exchange/transport rxns with no GPRs
%           *count(4):   #other rxns
% 
% Benjamin Sanchez. Last edited: 2017-03-05
% Ivan Domenzain.   Last edited: 2018-09-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model_data = getEnzymeCodes(model,action)

if nargin<2
    action = 'display';
end

%Standardize grRules to avoid wrong enzyme codes assignments to reactions
[grRules,~]   = standardizeGrRules(model,true);
model.grRules = grRules;

data      = load('../../databases/ProtDatabase.mat');
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

[m,n]      = size(model.S);
substrates = cell(n,20);
products   = cell(n,20);
uniprots   = cell(n,20);
Genes      = cell(n,20);
EC_numbers = cell(n,20);
MWs        = zeros(n,20);
isrev      = zeros(n,1);
count      = zeros(4,1);
rgmat      = full(model.rxnGeneMat);
conflicts  = cell(1,4);
for i = 1:n
    ks  = 1;
    kp  = 1;
    dir = 0;
    inv = 0;
    DB = '';
    %Save the substrates and products (if rxn is reversible):
    for j = 1:m
        if model.S(j,i) < 0 && model.ub(i) > 0
            substrates{i,ks} = model.metNames{j};
            ks  = ks+1;
            dir = 1;
        elseif model.S(j,i) > 0 && model.lb(i) < 0
            products{i,kp} = model.metNames{j};
            kp  = kp+1;
            inv = 1;
        end
    end
    %isrev(i) = 0 if rxn is blocked, = 1 if non-reversible, and = 2 if
    %reversible:
    isrev(i) = dir + inv;
    if ~isempty(model.grRules{i})
        %Find match in Swissprot:
        [new_uni,new_EC,new_MW,newGene,multGenes] = findInDB(model.grRules{i},DBprotSwissprot,DBgenesSwissprot,DBecNumSwissprot,DBMWSwissprot);
        if ~isempty(union_string(new_EC))
            count(1) = count(1) + isrev(i);
            DBase    = 'swissprot';
            if ~isempty(multGenes{1})
                multGenes{3} = DBase;
            end
        else
            %Find match in KEGG:
            [new_uni,new_EC,new_MW,newGene,multGenes] = findInDB(model.grRules{i},DBprotKEGG,DBgenesKEGG,DBecNumKEGG,DBMWKEGG);
            if ~isempty(union_string(new_EC))
                count(2) = count(2) + isrev(i);
                DBase    = 'kegg';
                if ~isempty(multGenes{1})
                    multGenes{3} = DBase;
                end
            else
                %Check if rxn is an exchange/transport rxn with no GPRs:
                GPRs       = sum(rgmat(i,:));
                rxn_name   = lower(model.rxnNames{i});
                exchange   = contains(rxn_name,'exchange');
                uptake     = contains(rxn_name,'uptake');
                production = contains(rxn_name,'production');
                transport  = contains(rxn_name,'transport');
                if (exchange || uptake || production || transport) && GPRs == 0
                    count(3) = count(3) + isrev(i);
                else
                    count(4) = count(4) + isrev(i);
                end
            end
        end
    
        for j = 1:length(new_uni)
            uniprots{i,j} = new_uni{j};
            Genes{i,j}    = newGene{j};
            if isempty(new_EC{j})
                EC_numbers{i,j} = union_string(new_EC);
            else
                EC_numbers{i,j} = new_EC{j};
            end
            MWs(i,j) = new_MW(j);
        end
        
        if ~isempty(multGenes{1})
            %Rxn index
            conflicts{1} = [conflicts{1};i];
            %Gene IDs
            conflicts{2} = [conflicts{2};multGenes{1}];
            %Indexes in DB
            conflicts{3} = [conflicts{3};multGenes{2}];
            %DB name
            conflicts{4} = [conflicts{4};multGenes{3}];
            if strcmpi(action,'add')
                if strcmpi(DBase,'swissprot')
                    [uni,EC,MW,Genes] = addMultipleMatches(uni,EC,MW,Genes,multGenes,swissprot);
                elseif strcmpi(DBase,'KEGG')
                    [uni,EC,MW,Genes] = addMultipleMatches(uni,EC,MW,Genes,multGenes,kegg);
                end
            end
        end
    end
    if rem(i,100) == 0 | i == n
        disp(['Getting enzyme codes: Ready with rxn ' int2str(i)])
    end
end
model_data.model        = model;
model_data.substrates   = substrates;
model_data.products     = products;
model_data.matchedGenes = Genes;
model_data.uniprots     = uniprots;
model_data.EC_numbers   = EC_numbers;
model_data.MWs          = MWs;
model_data.count        = count;
%Display error message with the multiple gene-protein matches found
if strcmpi(action,'display') && ~isempty(conflicts{1})
    displayErrorMessage(conflicts,swissprot,kegg)
end

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