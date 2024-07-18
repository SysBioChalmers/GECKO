classdef KEY_CLASSNAME < ModelAdapter
    methods
        function obj = KEY_CLASSNAME()
            % Set initial values of the obj.params - they can be changed by the user
            
            % Directory where all model-specific files and scripts are kept.
            % Is assumed to follow the GECKO-defined folder structure.
            obj.params.path = fullfile('KEY_PATH', 'KEY_NAME');

			% Path to the conventional GEM that this ecModel will be based on.
			obj.params.convGEM = fullfile(obj.params.path,'models','yourModel.xml');

			% Average enzyme saturation factor
			obj.params.sigma = 0.5;

			% Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.5;

			% Fraction of enzymes in the model [g enzyme/g protein]
			obj.params.f = 0.5;
            
            % Growth rate the model should be able to reach when not
            % constraint by nutrient uptake (e.g. max growth rate) [1/h]
			obj.params.gR_exp = 0.41;

			% Provide your organism scientific name
			obj.params.org_name = 'genus species';
            
            % Taxonomic identifier for Complex Portal
            obj.params.complex.taxonomicID = [];

			% Provide your organism KEGG ID, selected at
			% https://www.genome.jp/kegg/catalog/org_list.html
			obj.params.kegg.ID = 'sce';
            % Field for KEGG gene identifier; should match the gene
            % identifiers used in the model. With 'kegg', it takes the
            % default KEGG Entry identifier (for example YER023W here:
            % https://www.genome.jp/dbget-bin/www_bget?sce:YER023W).
            % Alternatively, gene identifiers from the "Other DBs" section
            % of the KEGG page can be selected. For example "NCBI-GeneID",
            % "UniProt", or "Ensembl". Not all DB entries are available for
            % all organisms and/or genes.
            obj.params.kegg.geneID = 'kegg';

			% Provide what identifier should be used to query UniProt.
            % Select proteome IDs at https://www.uniprot.org/proteomes/
            % or taxonomy IDs at https://www.uniprot.org/taxonomy.
            obj.params.uniprot.type = 'taxonomy'; % 'proteome' or 'taxonomy'
			obj.params.uniprot.ID = '559292'; % should match the ID type
            % Field for Uniprot gene ID - should match the gene ids used in the 
            % model. It should be one of the "Returned Field" entries under
            % "Names & Taxonomy" at this page: https://www.uniprot.org/help/return_fields
            obj.params.uniprot.geneIDfield = 'gene_oln';
            % Whether only reviewed data from UniProt should be considered.
            % Reviewed data has highest confidence, but coverage might be (very)
            % low for non-model organisms
            obj.params.uniprot.reviewed = false;

			% Reaction ID for glucose exchange reaction (or other preferred carbon source)
			obj.params.c_source = 'r_1714'; 

			% Reaction ID for biomass pseudoreaction
			obj.params.bioRxn = 'r_4041';

			% Name of the compartment where the protein pseudometabolites
            % should be located (all be located in the same compartment,
            % this does not interfere with them catalyzing reactions in
            % different compartments). Typically, cytoplasm is chosen.
			obj.params.enzyme_comp = 'cytoplasm';
        end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
            % Indicates how spontaneous reactions are identified. Here it
            % is done by the reaction have 'spontaneous' in its name.
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end