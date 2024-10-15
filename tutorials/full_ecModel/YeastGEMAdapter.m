classdef YeastGEMAdapter < ModelAdapter 
	methods
		function obj = YeastGEMAdapter()
			obj.params.path = fullfile(findGECKOroot,'tutorials','full_ecModel');

			obj.params.convGEM = fullfile(obj.params.path,'models','yeast-GEM.yml');

			obj.params.sigma = 0.5;

			obj.params.Ptot = 0.5;

			obj.params.f = 0.5;
			
			obj.params.gR_exp = 0.41;

			obj.params.org_name = 'saccharomyces cerevisiae';
			
			obj.params.complex.taxonomicID = 559292;

			obj.params.kegg.ID = 'sce';

			obj.params.kegg.geneID = 'kegg';

			obj.params.uniprot.type = 'proteome';

			obj.params.uniprot.ID = 'UP000002311';

			obj.params.uniprot.geneIDfield = 'gene_oln';

			obj.params.uniprot.reviewed = true;

			obj.params.c_source = 'r_1714'; 

			obj.params.bioRxn = 'r_4041';

            % Name of the compartment where the protein pseudometabolites
            % should be located (all be located in the same compartment,
            % this does not interfere with them catalyzing reactions in
            % different compartments). Typically, cytoplasm is chosen.
			obj.params.enzyme_comp = 'cytoplasm';			
		end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
		
	end
end
