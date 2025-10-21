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

            % Hyperparameters for SMC-ABC Bayesian kcat fitting, as
            % described in full_tutorial/protocol.m. Note: these are not
            % described in the Nature Protocols paper, but only introduced
            % from GECKO 3.3.0.
            % Number of sampled models per generation
            obj.params.bayesian.samplesPerGen       = 100;
            % Number of sampled models for the first generation. A larger
            % number is appropriate, as it allows to explore a wider range
            % of kcat values in the beginning
            obj.params.bayesian.samplesFirstGen     = 200;
            % The number of best performing models to keep during each
            % generation, to define the next prior distribution of kcat
            % values
            obj.params.bayesian.bestSamplesToKeep   = 80;
            % RMSE threshold when SMC-ABC should halt
            obj.params.bayesian.rmseThreshold       = 0.2;
            % Maximum number of generations after which SMC-ABC should halt
            obj.params.bayesian.maxGenerations      = 50;
		end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end
