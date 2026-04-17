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

            %% Hyperparameters for Bayesian kcat fitting
            % Default initial uncertainty (standard deviation in log-space) for kcat values
            obj.params.bayesian.sigma0logDefault    = 0.5;
            % Data sources for kcat values, ordered from least to most trusted
            obj.params.bayesian.kcatSources         = {'dlkcat','brenda','custom'};
            % Initial uncertainty for each source (lower = more trusted data)
            obj.params.bayesian.sigma0logSource     = [0.4; 0.2; 0.1];

            % Default shrinkage threshold: standard deviations required for full posterior update
            obj.params.bayesian.shrinkThrDefault    = 1.5;
            % Source-specific shrinkage thresholds (higher = more resistant to change)
            obj.params.bayesian.shrinkThrSource     = [1.5, 3.5, 5.5];
            % Default maximum posterior/prior variance ratio (prevents runaway uncertainty)
            obj.params.bayesian.varianceCapDefault  = 10;
            % Source-specific variance caps (tighter for more trusted sources)
            obj.params.bayesian.varianceCapSource   = [10,4,2];

            % Default threshold below which kcat is locked to prior (-1 = never lock)
            obj.params.bayesian.forcePriorThrDefault = -1;
            % Source-specific thresholds for locking to prior (higher = easier to lock)
            obj.params.bayesian.forcePriorThrSource  = [-1, 4, 8];
            % Minimum deviation (in σ units) required for parameter update (promotes sparsity)
            obj.params.bayesian.sparsityThreshold   = 0.3;

            % Generations at which to adjust sampling strategy
            obj.params.bayesian.scheduleGenerations = [1, 2, 9, 15];
            % Number of samples to draw at each scheduled generation (decreases over time)
            obj.params.bayesian.scheduleSamples     = [1000, 800, 600, 400];

            % RMSE percentile threshold for accepting samples (lower = stricter selection)
            obj.params.bayesian.targetAccept        = 10;
            % Minimum fraction of samples to retain each generation
            obj.params.bayesian.minKeep             = 0.3;
            % Maximum fraction of samples to retain each generation
            obj.params.bayesian.maxKeep             = 0.6;
            
            % Stop optimization when RMSE falls below this threshold
            obj.params.bayesian.rmseThreshold       = 0.2;
            % Maximum number of ABC-SMC generations before termination
            obj.params.bayesian.maxGenerations      = 150;
            % Post-optimization pruning
            obj.params.bayesian.enablePruning       = true; % Enable post-hoc sensitivity analysis
            obj.params.bayesian.prunRMSEtol         = 0.02; % Max acceptable RMSE increase (e.g., 2%)
            obj.params.bayesian.prunMinDev          = 1.0;  % Only test params that moved >1σ from prior
   end

        function ecModel = makeModelAnaerobic(obj,ecModel)
            % Taken from yeast-GEM 9.0.2
            ecModel = anaerobicModel_GECKO(ecModel);
        end
	    function ecModel = changeProteinBiomass(obj,ecModel,Ptot)
            % Skip this step for now
            ecModel = ecModel;
            % Taken from yeast-GEM 9.0.2
            %ecModel = scaleBioMass_GECKO(ecModel,'protein',Ptot,'carbohydrate',false);
        end
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end
