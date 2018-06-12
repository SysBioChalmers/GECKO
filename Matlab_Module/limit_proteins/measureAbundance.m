%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,count] = measureAbundance(model,abundance_file)
% 
%
% Benjamín J. Sánchez. Last edited: 2016-03-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,count] = measureAbundance(enzymes,abundance_file)
current = pwd;
%Read downloaded data of abundance:
fID       = fopen(abundance_file);
data      = textscan(fID,'%s %s %f','delimiter','\t','HeaderLines',12);
genes     = data{2};
abundance = data{3};
fclose(fID);

%Load KEGG data:
cd ../../Databases
data      = load('ProtDatabase.mat');
swissprot = data.swissprot;

%Main loop:
MW_ave  = mean(cell2mat(swissprot(:,5)));
Pmodel  = 0;
Ptot    = 0;
counter = zeros(size(genes));
for i = 1:length(genes)
    MW = MW_ave;
    %Trim gene name:
    gene_name = genes{i};
    gene_name = gene_name(strfind(gene_name,'.')+1:end);
    %Find gene in swissprot database:
    for j = 1:length(swissprot)
        if ~isempty(strfind(swissprot{j,3},gene_name))
            MW = swissprot{j,5};
            %Check if uniprot is in model:
            for k = 1:length(enzymes)
                if strcmp(enzymes{k},swissprot{j,1})
                    counter(i) = 1;
                    Pmodel = Pmodel + MW*abundance(i)/1000000;
                end
            end
        end
    end
    Ptot = Ptot + MW*abundance(i)/1000000;
    disp(['Calculating total abundance: Ready with gene ' gene_name])
end

f     = Pmodel/Ptot;
count = [length(counter);sum(counter(:,1))];
cd (current)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%