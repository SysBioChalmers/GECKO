classdef HumanGEMAdapter < ModelAdapter 
    methods
        function obj = HumanGEMAdapter()
            %Set initial values of the parameters - they can be changed by the user
            
            %these parameters are just copied from getParams in ecModels
			
            %Average enzyme saturation factor
			obj.params.sigma = 0.5;

			%Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.505717472;  %Average across NCI60 cell lines

			%Minimum growth rate the model should grow at [1/h]
			obj.params.gR_exp = 0.020663429; %[g/gDw h]/Average across NCI60 cell lines

			%Provide your organism scientific name
			obj.params.org_name = 'homo sapiens';

			%Provide your organism KEGG ID
			obj.params.keggID = 'hsa';
            
            %Taxon id for Uniprot
            obj.params.taxonID = '9606';
            
            %Field for Uniprot gene id - should match the gene ids used in the 
            %GPRs. Note that this is a field in the web request to uniprot - 
            %it has to match one of the fields there
            %So, Ensembl gives ensembl transcripts, not genes, so we need to
            %convert the GPRs and match them to gene ids
            obj.params.uniprotGeneIdField = 'gene_primary';

			%The name of the exchange reaction that supplies the model with carbon (rxnNames)
			obj.params.c_source = 'MAR09034'; 

			%Experimental carbon source uptake (optional)
			obj.params.c_UptakeExp = 0.641339301; %[mmol/gDw h]/Average across NCI60 cell lines

			%Rxn Id for non-growth associated maitenance pseudoreaction
			obj.params.NGAM = 'MAR09931'; %this is not used

			%Compartment name in which the added enzymes should be located
			obj.params.enzyme_comp = 'Cytosol';
        end
		
        function result = getFilePath(obj, filename)
			result = filename; % TODO: Look at how this should be solved - look at the GeckoLight solution below
			%result = strcat(GeckoLightInstall.getGeckoLightMainPath(), 'data/humanGEM/', filename);
		end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			rxns_tsv = importTsvFile(strcat(getHumanGEMRootPath(),'model/reactions.tsv'));
			spont = rxns_tsv.spontaneous;
			spontRxnNames = rxns_tsv.rxns;
        end
		
        %Ensembl gene ids are not available in uniprot (ensembl returns transcripts, not genes)
        %Therefore, we collect gene symbols from uniprot, and need to convert the ensembl genes
        %here to gene symbols as well
        function genes = getUniprotCompatibleGenes(obj,model)
            % get the path
            tmpfile = fullfile(getHumanGEMRootPath(),'model','genes.tsv');

            % import as structure, convert to table, and extract header
            tmp = struct2table(importTsvFile(tmpfile));
            conv_key = table2array(tmp);
            clear tmp

            %convert genes
            genes = model.genes;
            [~,ia,ib] = intersect(genes, conv_key(:,1));
            genes(ia) = conv_key(ib,5);
            
            %A special canse - GGT2 is not found in uniprot - however GGT2P is
            %This is supposed to be a pseudogene, but I suspect the sequence is very similar, so let's use that
            genes{strcmp(genes, 'GGT2')} = 'GGT2P';
            %Also replace empty strings with something else to avoid matches on empty strings in uniprot
            genes(strcmp(genes, '')) = {'<<<EMPTY>>>'};
        end

        
		function model = manualModifications(obj,model) %default is to do nothing
			%So, there are some of these in ecModels - it is a bit unclear if any of these are relevant here
			%we do nothing for now.
			%In general, manual modifications should be done to the model before sending it in.
		end
	end
end
