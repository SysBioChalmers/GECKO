classdef YeastGEMAdapter < ModelAdapter 
	methods
		function obj = YeastGEMAdapter()
			obj.params.path = fullfile(findGECKOroot,'tutorials','full_ecModel');
            addpath(fullfile(obj.params.path,'code'));

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
            obj.params.bayesian.initSDmultiplDef    = 1.5; % Default multiplier to define kcat standard deviation (kcat * initSDmultiplDef)
            obj.params.bayesian.kcatSources         = {'brenda','dlkcat'}; % kcat sources for which alternative multipliers are defined
            obj.params.bayesian.initSDmultipl       = [0.15; 0.6]; % Multipliers that overwrite initSDmultplDef, matching kcatSources
            obj.params.bayesian.samples1            = 1500; % Number of kcat samples in generation 1
            obj.params.bayesian.samples2_5          = 500; % Number of kcat samples in generation 2-5
            obj.params.bayesian.samples6_end        = 250; % Number of kcat samples in generation 6 onwards
            obj.params.bayesian.targetAccept1       = 45; % RMSE percentile to accept in generation 1
            obj.params.bayesian.targetAccept2_5     = 35; % RMSE percentile to accept in generation 2-5
            obj.params.bayesian.targetAccept6_end   = 25; % RMSE percentile to accept in generation 6 onwards
            obj.params.bayesian.minKeep             = 0.05; % Minimum fraction of samples to keep
            obj.params.bayesian.maxKeep             = 0.30; % Maximum fraction of samples to keep
            obj.params.bayesian.rMax                = 200; % Hard cap on PCA rank (20–300 typical)
            obj.params.bayesian.cExpl               = 1.25; % inflate subspace std for exploration (1.1–1.5)
            obj.params.bayesian.tauResidual         = 1e-3; % residual (isotropic) variance in log-space (1e-4–1e-2
            obj.params.bayesian.rmseThreshold       = 0.2; % RMSE threshold to halt and output best posterior kcats
            obj.params.bayesian.maxGenerations      = 200; % Maximum number of generations before returning best posterior kcats
        end

        function ecModel = makeModelAnaerobic(ecModel)
            % Taken from yeast-GEM 9.0.2
            ecModel = anaerobicModel_GECKO(ecModel);
        end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end
