classdef HumanGEMAdapter < ModelAdapter 
    methods
        function obj = HumanGEMAdapter()
            geckoPath = findGECKOroot;
            obj.params.path = fullfile(geckoPath,'tutorials','ecHumanGEM');

			obj.params.convGEM = fullfile(HumanGEMAdapter.getHumanGEMRootPath(),'model','Human-GEM.xml');

			obj.params.sigma = 0.1; %This was changed to a low number to give a reasonable growth rate - this should be investigated more

			obj.params.Ptot = 0.505717472;  %Average across NCI60 cell lines

			obj.params.f = 0.412; %estimated as TPM of model genes / all genes in GTEx

			obj.params.gR_exp = 0.020663429; %[g/gDw h]/Average across NCI60 cell lines
            
			obj.params.c_UptakeExp = 0.641339301; %[mmol/gDw h]/Average across NCI60 cell lines
            
			obj.params.org_name = 'homo sapiens';
            
            obj.params.complex_org_name = 'Homo sapiens';

			obj.params.keggID = 'hsa';

            obj.params.keggGeneIdentifier = 'Ensembl';

            obj.params.uniprotIDtype = 'proteome';

			obj.params.uniprotID = 'UP000005640';

            obj.params.uniprotGeneIdField = 'gene_primary';

            obj.params.uniprotReviewed = true;
            
			obj.params.c_source = 'MAR09034'; 

			obj.params.bioRxn = 'MAR13082';

			obj.params.NGAM = 'MAR09931'; %This is not used

			obj.params.enzyme_comp = 'Cytosol';
            
        end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			rxns_tsv = importTsvFile(strcat(HumanGEMAdapter.getHumanGEMRootPath(),'model/reactions.tsv'));
			spont = rxns_tsv.spontaneous;
			spontRxnNames = rxns_tsv.rxns;
        end
		
        % Ensembl gene ids are not available in uniprot (ensembl returns transcripts, not genes)
        % Therefore, we collect gene symbols from uniprot, and need to convert the ensembl genes
        % here to gene symbols as well
        function genes = getUniprotCompatibleGenes(obj,inGenes)
            % Get the path
            tmpfile = fullfile(HumanGEMAdapter.getHumanGEMRootPath(),'model','genes.tsv');

            % Import as structure, convert to table, and extract header
            tmp = struct2table(importTsvFile(tmpfile));
            conv_key = table2array(tmp);
            clear tmp

            % Convert genes
            genes = inGenes;
            [~,ia,ib] = intersect(genes, conv_key(:,1));
            genes(ia) = conv_key(ib,5);
            
            % A special canse - GGT2 is not found in uniprot - however GGT2P is
            % This is supposed to be a pseudogene, but I suspect the sequence is very similar, so let's use that
            genes{strcmp(genes, 'GGT2')} = 'GGT2P';
            % Also replace empty strings with something else to avoid matches on empty strings in uniprot
            genes(strcmp(genes, '')) = {'<<<EMPTY>>>'};
        end
    end
    
    methods(Static)
        function path = getHumanGEMRootPath()
            path = fileparts(which('Human-GEM.mat'));% This file should be on the path
            path = strrep(path, '\', '/'); % get rid of backslashes in Windows
            if ~endsWith(path, '/')
                path = strcat(path,'/');
            end
            % Now remove the model/ at the end
            path = (path(1:strlength(path)-6));
        end
	end
end
