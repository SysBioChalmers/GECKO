function [uni,EC,MW,Genes,conflicts] = findInDB(grRule,DBprot,DBgenes,DBecNum,DBMW, geneIndex, geneHashMap)
% findInDB
%   Gets the uniprots and EC numbers for a given rxn into a given database
%
%   grRule     Genes association for a given metabolic reaction
%   DBprot     array of uniprot IDs, taken from swissprot or kegg database
%   DBgenes    array of genes corresponding to DBprot
%   DBecNum    array of EC numbers corresponding to DBprot
%   DBMW       array of molecular weights in Da corresponding to DBprot
%
%   uni        Array containing all the uniprot IDs related to an isoenzyme,
%              for a given reaction in each of its cells.
%   EC         Array containing all the EC numbers related to an isoenzyme,
%              for a given reaction in each of its cells.
%   Genes      Array containing all the genes related to an isoenzyme,
%              for a given reaction in each of its cells.
%   conflicts  Contains those genes with multiple protein matches. 
%
%   Usage: [uni,EC,MW,Genes,conflicts] = findInDB(grRule,DBprot,DBgenes,DBecNum,DBMW)
%
%   Benjamin J. Sanchez, 2017-08-10
%   Ivan Domenzain,      2018-09-06
%	Johan Gustafsson, 	 2021-07-02 introduced optimizations from GeckoLight
%

%Find simple gene sets for reaction:
gene_sets = getSimpleGeneSets(grRule);
%Preallocate function variables
uni       = cell(size(gene_sets));
EC        = cell(size(gene_sets));
MW        = zeros(size(gene_sets));
Genes     = cell(size(gene_sets));
conflicts = cell(1,2);

for i = 1:length(gene_sets)
    %Split the gene set and match each gene:
    gene_set  = strsplit(gene_sets{i},' and '); %this string splitting together with getSimpleGeneSets takes 2.6 s, that is ok.
    uni_set   = cell(size(gene_set));
    uni_set_idx   = cell(size(gene_set));
    EC_set    = cell(size(gene_set));
    genesCell = cell(size(gene_set));
    
    %find index indices of all genes in one call

    %this loop is slow
    for j = 1:length(gene_set)
         gene = gene_set{j};
         matches = [];
         %Get the indexes for all of the proteins related to gene
         if isKey(geneHashMap,gene) %annoyingly, this seems to be needed
             matchInd = cell2mat(values(geneHashMap,gene_set(j)));
             matches = geneIndex{matchInd};
             %matches2=find(sum(strcmpi(gene,DBgenes),2));
             %if (length(matches) ~= length(matches2))
             %   disp('misfit length')
             %end
             %for i = 1:length(matches)
             %   if (~all(matches == matches2))
             %       disp('misfit numbers')
             %   end
             %end
         end
%    end %remove later
         %matches=find(sum(strcmpi(gene,DBgenes),2));
         if ~isempty(matches)
            uni_set{j}   = DBprot{matches};
            uni_set_idx{j} = matches;
            genesCell{j} = gene;
            %Get the indexes of the matched protein(s) with non-empty EC#s
           nonEmpty = matches(~cellfun(@isempty,DBecNum(matches)));
            if ~isempty(nonEmpty)
                [geneECs,ia] = unique(DBecNum(nonEmpty),'stable');
                %For genes with multiple proteins associated to several
                %non-empty ecNumbers
                if length(geneECs)>1
                    indexes = nonEmpty(ia);
                    ecNum = DBecNum(indexes);
                    %Save first match
                    uni_set{j} = DBprot{indexes(1)};
                    uni_set_idx{j} = indexes(1);
                    EC_set{j}  = getECstring(EC_set{j},ecNum{1});
                    %Save additional matches as potential conflicts
                    if ~ismember(gene,conflicts{1})
                        conflicts{1} = [conflicts{1}; {gene}];
                        conflicts{2} = [conflicts{2}; {indexes}]; 
                    end
                else
                    %If there is a single unique ec number for the matched
                    %protein(s), then choose the lightest protein
                    [~,minW]   = min(cell2mat(DBMW(nonEmpty)));
                    matches    = nonEmpty(minW);
                    uni_set{j} = DBprot{matches};
                    uni_set_idx{j} = matches;
                    EC_set{j}  = getECstring(EC_set{j},DBecNum{matches});
                end
            end
        else
            uni_set{j} = '';
            uni_set_idx{j} = NaN;
        end
                             
        if isempty(EC_set{j})
            EC_set{j} = '';
        end
     end
    %Uniprot: Delete repeated and empty spaces
    [uni_set,repeated] = deleteRepeated(uni_set);
    uni_set_idx        = uni_set_idx(~repeated);
    genesCell          = genesCell(~repeated);
    emptyIndexes       = cellfun(@isempty,uni_set);
    uni_set            = uni_set(~emptyIndexes);
    uni_set_idx        = uni_set_idx(~emptyIndexes);
    genesCell          = genesCell(~emptyIndexes);

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
        %index2     = strcmpi(DBprot(:),uni_set{j});
        index = uni_set_idx{j};
        %if (length(find(index2)) ~= length(index))
        %    disp('Length error')
        %elseif (~all(index == find(index2)))
        %    disp('value error')
        %end
        minWeight = min(DBMW{index});
        MW(i)     = MW(i) + minWeight;
    end
    %Add new codes as new possible isoenzymes:
    uni{i}   = strjoin(uni_set,' ');
    EC{i}    = strjoin(EC_set,' ');
    Genes{i} = strjoin(genesCell,' ');
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