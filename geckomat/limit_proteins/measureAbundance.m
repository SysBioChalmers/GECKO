%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,count] = measureAbundance(enzymes)
% 
% Benjamin J. Sanchez. Last edited: 2018-08-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,count] = measureAbundance(enzymes)

%Read downloaded data of abundance:
fID       = fopen('../../databases/prot_abundance.txt');
data      = textscan(fID,'%s %s %f','delimiter','\t','HeaderLines',12);
genes     = data{2};
genes     = strrep(genes,'4932.','');
abundance = data{3};
fclose(fID);

%Load swissprot data:
data      = load('../../databases/ProtDatabase.mat');
swissprot = data.swissprot;
for i = 1:length(swissprot)
    swissprot{i,3} = strsplit(swissprot{i,3},' ');
end

%Main loop:
MW_ave  = mean(cell2mat(swissprot(:,5)));
concs   = zeros(size(genes));
counter = false(size(genes));
for i = 1:length(genes)
    MW = MW_ave;
    %Find gene in swissprot database:
    for j = 1:length(swissprot)
        if sum(strcmp(swissprot{j,3},genes{i})) > 0
            MW = swissprot{j,5};	%g/mol
            %Check if uniprot is in model:
            if sum(strcmp(enzymes,swissprot{j,1})) > 0
                counter(i) = true;
            end
        end
    end
    concs(i) = MW*abundance(i);     %g/mol(tot prot)
    if rem(i,100) == 0
        disp(['Calculating total abundance: Ready with ' num2str(i) '/' ...
              num2str(length(genes)) ' genes '])
    end
end

f     = sum(concs(counter))/sum(concs);
count = [length(counter);sum(counter)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
