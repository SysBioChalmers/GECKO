classdef YeastGEMAdapter < ModelAdapter 
	methods
		function obj = YeastGEMAdapter()
			obj.params.path = fullfile(findGECKOroot,'tutorials','tutorial_yeast-GEM');

			obj.params.convGEM = fullfile(obj.params.path,'models','yeast-GEM.xml');

			obj.params.sigma = 0.5;

			obj.params.Ptot = 0.5;

			obj.params.f = 0.5;
			
			obj.params.gR_exp = 0.41;

			obj.params.org_name = 'saccharomyces cerevisiae';
			
			obj.params.complex_org_name = 'Saccharomyces cerevisiae';

			obj.params.keggID = 'sce';

			obj.params.keggGeneIdentifier = 'kegg';

			obj.params.uniprotIDtype = 'taxonomy_id';

			obj.params.uniprotID = '559292';

			obj.params.uniprotGeneIdField = 'gene_oln';

			obj.params.uniprotReviewed = true;

			obj.params.c_source = 'r_1714'; 

			obj.params.bioRxn = 'r_4041';

			obj.params.NGAM = 'r_4046';

			obj.params.enzyme_comp = 'cytoplasm';

			obj.params.exch_names{1} = 'growth';
			obj.params.exch_names{2} = 'D-glucose exchange (reversible)';
			obj.params.exch_names{3} = 'oxygen exchange (reversible)';
			obj.params.exch_names{4} = 'carbon dioxide exchange';

			obj.params.bio_comp{1} = 'protein';
			obj.params.bio_comp{2} = 'carbohydrate';
			obj.params.bio_comp{3} = 'lipid backbone';
			obj.params.bio_comp{4} = 'lipid chain';

			obj.params.pol_cost(1) = 37.7; % Ptot 
			obj.params.pol_cost(2) = 12.8; % Ctot
			obj.params.pol_cost(3) = 26.0; % RNA 
			obj.params.pol_cost(4) = 26.0; % DNA

			obj.params.oxPhos{1} = 'r_1021';
			obj.params.oxPhos{2} = 'r_0439';
			obj.params.oxPhos{3} = 'r_0438';
			obj.params.oxPhos{4} = 'r_0226';
			
		end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
		
	end
end
