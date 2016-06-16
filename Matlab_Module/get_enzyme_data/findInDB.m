%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [uni,EC] = findInDB(rxn_pos,model,DB)
% Matches the uniprot and EC number for a given rxn into a given database.
%
% Benjamín J. Sánchez. Last edited: 2015-08-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uni,EC,MW] = findInDB(rxn_pos,model,DB)

%Find simple gene sets for reaction:
[gene_sets,~] = getAllPath(model,model.rxns(rxn_pos));
uni           = cell(size(gene_sets));
EC            = cell(size(gene_sets));
MW            = zeros(size(gene_sets));

for i = 1:length(gene_sets)
    %Split the gene set and match each gene:
    gene_set = strsplit(gene_sets{i},' AND ');
    uni_set  = cell(size(gene_set));
    EC_set   = cell(size(gene_set));
    for j = 1:length(gene_set)
        for k = 1:length(DB)
            if ~isempty(strfind(DB{k,3},gene_set{j}))
                uni_set{j} = DB{k,1};
                if ~isempty(DB{k,4})
                    %Always prefer protein with EC value:
                    uni_set{j} = DB{k,1};
                    new_EC_set = strsplit(DB{k,4},' ');
                    for l = 1:length(new_EC_set)
                        EC_set{j}  = ['EC' new_EC_set{l} ' '];
                    end
                end
            end
        end
        if isempty(EC_set{j})
            EC_set{j} = '';
        else
            EC_set{j} = EC_set{j}(1:end-1);
        end
    end
    %Uniprot: Delete repeated and empty spaces
    [uni_set,~] = deleteRepeated(uni_set);
    uni_set     = uni_set(~cellfun('isempty',uni_set));
    
    %EC: Find union and intersection between all units (only applies for
    %complexes, i.e. length(EC_set) > 1):
    uni_EC = strsplit(EC_set{1},' ');
    int_EC = uni_EC;
    for j = 2:length(EC_set)
        other_EC = strsplit(EC_set{j},' ');
        if isempty(uni_EC)
            uni_EC = other_EC;
            int_EC = other_EC;
        elseif ~isempty(other_EC)
            uni_EC = compare_wild([uni_EC other_EC]);
            int_EC = intersection(int_EC,other_EC);
        end
    end
    %Use the intersection found, if any. If not, use the union: 
    if isempty(int_EC)
        EC_set = uni_EC;
    else
        EC_set = int_EC;
    end
    
    %MW: use uniprot codes + swissprot table:
    for j = 1:length(uni_set)
        for k = 1:length(DB)
            if strcmpi(uni_set{j},DB{k,1})
                MW(i) = MW(i) + DB{k,5};
            end
        end
    end
    
    %Add new codes as new possible isoenzymes:
    uni{i} = strjoin(uni_set,' ');
    EC{i}  = strjoin(EC_set,' ');
end

%Delete repeated Uniprots (and the corresponding ECs):
[uni,deleted] = deleteRepeated(uni);
EC(deleted)   = [];
MW(deleted)   = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_EC = intersection(prev_EC,new_EC)
%Finds the common elements between two cell arrays, if any. Also considers
%wildcards (e.g. if 'EC1.1.1.1' is in one array and 'EC1.1.1.-' is in the
%other one, then 'EC1.1.1.1' is added to the intersection).

int_EC = {};
for i = 1:length(prev_EC)
    for j = 1:length(new_EC)
        new_int = compare_wild({prev_EC{i} new_EC{j}});
        if length(new_int) == 1
            int_EC = [int_EC new_int];
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_EC = compare_wild(EC)
%Goes through a cell array of EC numbers, and erases any repetitions,
%considering also wildcards (e.g. will erase 'EC1.2.3.-' if 'EC1.2.3.4' is
%already present).

%Trim all EC numbers of wild cards (e.g. 'EC1.2.-.-' -> 'EC1.2.'):
EC_trimmed = EC;
for i = 1:length(EC)
    %Add a final dot for avoiding issues (e.g. 'EC1.1.1.1' & 'EC1.1.1.12'):
    ECi = [EC{i} '.'];
    pos = strfind(ECi,'-');
    if ~isempty(pos)
        ECi = ECi(1:pos(1)-1);
    end
    EC_trimmed{i} = ECi;
end

%Compare all EC numbers between them to find if 1 fits in the other:
non_repeated = true(1,length(EC));
for i = 1:length(EC)-1
    for j = i+1:length(EC)
        ECi = EC_trimmed{i};
        ECj = EC_trimmed{j};
        %If ECj fits in ECi then ECj can be disregarded:
        if strfind(ECi,ECj)
            non_repeated(j) = false;
        %Else, if ECi fits in ECj then ECi can be disregarded:
        elseif strfind(ECj,ECi)
            non_repeated(i) = false;
        end
    end
end

int_EC = EC(non_repeated);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%