classdef YourAdapter < ModelAdapter
    methods
        function obj = YourAdapter()
            % Set initial values of the parameters - they can be changed by the user
            
            % Directory where all model-specific files and scripts are kept.
            % Is assumed to follow the GECKO-defined folder structure. The
            % code below refers to tutorials/ecHumanGEM in the GECKO path.
            geckoPath = findGECKOroot;
            obj.params.path = fullfile(geckoPath,'tutorials','blankExample');

			% Path to the conventional GEM that this ecModel will be based on.
			obj.params.convGEM = fullfile(HumanGEMAdapter.getHumanGEMRootPath(),'model','yourModel.xml');

            % Average enzyme saturation factor
			obj.params.sigma = 0.1; %This was changed to a low number to give a reasonable growth rate - this should be investigated more

			% Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.505717472;  %Average across NCI60 cell lines

			% Fraction of enzymes in the model [g enzyme/g protein]
			obj.params.f = 0.412; %estimated as TPM of model genes / all genes in GTEx

			% Minimum growth rate the model should grow at [1/h]
			obj.params.gR_exp = 0.020663429; %[g/gDw h]/Average across NCI60 cell lines
            
			% Experimental carbon source uptake (optional)
			obj.params.c_UptakeExp = 0.641339301; %[mmol/gDw h]/Average across NCI60 cell lines
            
            % Provide your organism scientific name
			obj.params.org_name = 'homo sapiens';
            
            % Matching name for Complex Portal
            obj.params.complex_org_name = 'Homo sapiens';

			% Provide your organism KEGG ID, selected at
			% https://www.genome.jp/kegg/catalog/org_list.html
			obj.params.keggID = 'hsa';
            % Field for KEGG gene identifier; should match the gene
            % identifiers used in the model. With 'kegg', it takes the
            % default KEGG Entry identifier (for example YER023W here:
            % https://www.genome.jp/dbget-bin/www_bget?sce:YER023W).
            % Alternatively, gene identifiers from the "Other DBs" section
            % of the KEGG page can be selected. For example "NCBI-GeneID",
            % "UniProt", or "Ensembl". Not all DB entries are available for
            % all organisms and/or genes.
            obj.params.keggGeneIdentifier = 'Ensembl';

			% Provide what identifier should be used to query UniProt.
            % Select proteome IDs at https://www.uniprot.org/proteomes/
            % or taxonomy IDs at https://www.uniprot.org/taxonomy.
            obj.params.uniprotIDtype = 'proteome'; % 'proteome' or 'taxonomy_id'
			obj.params.uniprotID = 'UP000005640'; % should match the ID type
            
            % Field for Uniprot gene id - should match the gene ids used in the 
            % model. It should be one of the "Returned Field" entries under
            % "Names & Taxonomy" at this page: https://www.uniprot.org/help/return_fields
            % Human-GEM uses Ensemble IDs, but these are not easily retrievable
            % from Uniprot. Instead, we will gather gene_primary entries, and in
            % a separate function (getUniprotCompatibleGenes) are these matched
            % to Ensembl IDs.
            obj.params.uniprotGeneIdField = 'gene_primary';

            % Whether only reviewed data from UniProt should be considered.
            % Reviewed data has highest confidence, but coverage might be (very)
            % low for non-model organisms
            obj.params.uniprotReviewed = true;
            
			% Reaction ID for glucose exchange reaction (or other preferred carbon source)
			obj.params.c_source = 'SET_C_SOURCE'; 

			% Reaction ID for biomass pseudoreaction
			obj.params.bioRxn = 'SET_BIORXN';

			% Reaction ID for non-growth associated maitenance pseudoreaction
			obj.params.NGAM = 'SET_NGAM'; %This is not used

			% Compartment name in which the added enzymes should be located
			obj.params.enzyme_comp = 'SET_ENZYME_COMP';
            
        end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			rxns_tsv = importTsvFile(strcat(HumanGEMAdapter.getHumanGEMRootPath(),'model/reactions.tsv'));
			spont = rxns_tsv.spontaneous;
			spontRxnNames = rxns_tsv.rxns;
        end
		
    end
    
end
