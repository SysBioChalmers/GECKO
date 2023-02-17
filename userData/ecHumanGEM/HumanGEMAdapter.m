classdef HumanGEMAdapter < ModelAdapter 
    methods
        function obj = HumanGEMAdapter()
            %Set initial values of the parameters - they can be changed by the user
            
            %Directory where all model-specific files and scripts are kept.
            %Is assumed to follow the GECKO-defined folder structure. The
            %code below refers to userData/ecYeastGEM in the GECKO path.
            geckoPath = findGECKOroot;
            obj.params.path = fullfile(geckoPath,'userData','ecHumanGEM');

			%Path to the conventional GEM that this ecModel will be based on.
			obj.params.convGEM = fullfile(obj.params.path,'models','Human-GEM.xml');
            
            %Average enzyme saturation factor
			obj.params.sigma = 0.5;

			%Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.505717472;  %Average across NCI60 cell lines

			%Minimum growth rate the model should grow at [1/h]
			obj.params.gR_exp = 0.020663429; %[g/gDw h]/Average across NCI60 cell lines

			%Provide your organism scientific name
			obj.params.org_name = 'homo sapiens';
            
            %Matching name for Complex Portal
            obj.params.complex_org_name = 'Homo sapiens';

			%Provide your organism KEGG ID, selected at
			%https://www.genome.jp/kegg/catalog/org_list.html
			obj.params.keggID = 'hsa';

			%Provide your organism UniProt proteome, selected at
			%https://www.uniprot.org/proteomes/
			obj.params.uniprotID = 'UP000005640';
            
            %Field for Uniprot gene id - should match the gene ids used in the 
            %GPRs. Note that this is a field in the web request to uniprot - 
            %it has to match one of the fields there
            %So, Ensembl gives ensembl transcripts, not genes, so we need to
            %convert the GPRs and match them to gene ids
            obj.params.uniprotGeneIdField = 'gene_primary';

            %Whether only reviewed data from UniProt should be considered.
            %Reviewed data has highest confidence, but coverage might be (very)
            %low for non-model organisms
            obj.params.uniprotReviewed = true;
            
            %The name of the exchange reaction that supplies the model with carbon (rxnNames)
			obj.params.c_source = 'MAR09034'; 

            %The name of the exchange reaction that supplies the model with carbon (rxnNames)
			obj.params.bioRxn = 'MAR13082'; 

			%Experimental carbon source uptake (optional)
			obj.params.c_UptakeExp = 0.641339301; %[mmol/gDw h]/Average across NCI60 cell lines

			%Rxn Id for non-growth associated maitenance pseudoreaction
			obj.params.NGAM = 'MAR09931'; %this is not used

			%Compartment name in which the added enzymes should be located
			obj.params.enzyme_comp = 'Cytosol';
            
            %The pool size, fitted to NCI60 cell lines
            obj.params.standardProtPoolSize = 22.38315; %This is the value used in Gecko light
        end
		
        function result = getFilePath(obj, filename)
			result = filename; % TODO: Look at how this should be solved - look at the GeckoLight solution below
			%result = strcat(GeckoLightInstall.getGeckoLightMainPath(), 'data/humanGEM/', filename);
        end
        
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			rxns_tsv = importTsvFile(strcat(HumanGEMAdapter.getHumanGEMRootPath(),'model/reactions.tsv'));
			spont = rxns_tsv.spontaneous;
			spontRxnNames = rxns_tsv.rxns;
        end
		
        %Ensembl gene ids are not available in uniprot (ensembl returns transcripts, not genes)
        %Therefore, we collect gene symbols from uniprot, and need to convert the ensembl genes
        %here to gene symbols as well
        function genes = getUniprotCompatibleGenes(obj,inGenes)
            % get the path
            tmpfile = fullfile(HumanGEMAdapter.getHumanGEMRootPath(),'model','genes.tsv');

            % import as structure, convert to table, and extract header
            tmp = struct2table(importTsvFile(tmpfile));
            conv_key = table2array(tmp);
            clear tmp

            %convert genes
            genes = inGenes;
            [~,ia,ib] = intersect(genes, conv_key(:,1));
            genes(ia) = conv_key(ib,5);
            
            %A special canse - GGT2 is not found in uniprot - however GGT2P is
            %This is supposed to be a pseudogene, but I suspect the sequence is very similar, so let's use that
            genes{strcmp(genes, 'GGT2')} = 'GGT2P';
            %Also replace empty strings with something else to avoid matches on empty strings in uniprot
            genes(strcmp(genes, '')) = {'<<<EMPTY>>>'};
        end
    end
    
    methods(Static)
        function path = getHumanGEMRootPath()
            path = fileparts(which('Human-GEM.mat'));%This file should be on the path
            path = strrep(path, '\', '/'); %get rid of backslashes in Windows
            if ~endsWith(path, '/')
                path = strcat(path,'/');
            end
            %Now remove the model/ at the end
            path = (path(1:strlength(path)-6));
        end
	end
end
