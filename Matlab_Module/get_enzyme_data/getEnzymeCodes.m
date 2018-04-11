%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_data = getEnzymeCodes(model)
% Retrieves the enzyme codes for each of the reactions for a given genome
% scale model (GEM).
%
% INPUT:    a GEM file (.mat format)
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
% Ivan Domenzain.   Last edited: 2018-03-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model_data = getEnzymeCodes(model)

%Standardize grRules to avoid wrong enzyme codes assignments to reactions
[grRules,~]   = standardizeGrRules(model,true);
model.grRules = grRules;

cd ../../Databases
data      = load('ProtDatabase.mat');
swissprot = data.swissprot;
kegg      = data.kegg;
cd ../Matlab_Module/get_enzyme_data

swissprot = standardizeDatabase(swissprot);
kegg      = standardizeDatabase(kegg);

[m,n]      = size(model.S);
substrates = cell(n,20);
products   = cell(n,20);
uniprots   = cell(n,20);
EC_numbers = cell(n,20);
MWs        = zeros(n,20);
isrev      = zeros(n,1);
count      = zeros(4,1);
rgmat      = full(model.rxnGeneMat);

for i = 1:n
    ks  = 1;
    kp  = 1;
    dir = 0;
    inv = 0;
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
        [new_uni,new_EC,new_MW] = findInDB(model.grRules{i},swissprot);
        if ~isempty(union_string(new_EC))
            count(1) = count(1) + isrev(i);
        else
            %Find match in KEGG:
            [new_uni,new_EC,new_MW] = findInDB(model.grRules{i},kegg);
            if ~isempty(union_string(new_EC))
                count(2) = count(2) + isrev(i);
            else
                %Check if rxn is an exchange/transport rxn with no GPRs:
                GPRs       = sum(rgmat(i,:));
                rxn_name   = lower(model.rxnNames{i});
                exchange   = ~isempty(strfind(rxn_name,'exchange'));
                uptake     = ~isempty(strfind(rxn_name,'uptake'));
                production = ~isempty(strfind(rxn_name,'production'));
                transport  = ~isempty(strfind(rxn_name,'transport'));
                if (exchange || uptake || production || transport) && GPRs == 0
                    count(3) = count(3) + isrev(i);
                else
                    count(4) = count(4) + isrev(i);
                end
            end
        end
    
        for j = 1:length(new_uni)
            uniprots{i,j} = new_uni{j};
            if isempty(new_EC{j})
                EC_numbers{i,j} = union_string(new_EC);
            else
                EC_numbers{i,j} = new_EC{j};
            end
            MWs(i,j) = new_MW(j);
        end
    end
    disp(['Getting enzyme codes: Ready with rxn ' int2str(i)])
end

model_data.model      = model;
model_data.substrates = substrates;
model_data.products   = products;
model_data.uniprots   = uniprots;
model_data.EC_numbers = EC_numbers;
model_data.MWs        = MWs;
model_data.count      = count;

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