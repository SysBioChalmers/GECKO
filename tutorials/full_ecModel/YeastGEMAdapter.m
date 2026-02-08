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

            % Parameters for Bayesian kcat fitting
            obj.params.bayesian.samplesPerGen       = 150;
            obj.params.bayesian.samplesFirstGen     = 200;
            obj.params.bayesian.targetAccept        = 10; % RMSE percentile to accept in each iteration
            obj.params.bayesian.minKeep             = 0.3; % Minimum fraction of samples to keep
            obj.params.bayesian.maxKeep             = 0.6; % Maximum fraction of samples to keep
            obj.params.bayesian.alpha               = 0.7; % Exploit fraction
            obj.params.bayesian.cExpl               = 3.0; % Exploration inflation
            obj.params.bayesian.freezeStage         = 4; % Start freezing scale after this generation
            obj.params.bayesian.sigmaFloorFrac      = 0.10; % Keep ≥10% of initial multiplicative uncertainty
            obj.params.bayesian.rMax                = 150; % Cap PCA rank
            obj.params.bayesian.tauResidual         = 0.10; % Residual variance outside PCA subspace
            obj.params.bayesian.rmseThreshold       = 0.2;
            obj.params.bayesian.maxGenerations      = 200;
		end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end
