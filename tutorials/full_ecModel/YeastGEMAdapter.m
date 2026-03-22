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
            % Default initial stdev of the kcat log-normal distribution
            obj.params.bayesian.sigma0logDefault    = 0.5;                    
            % If other (smaller?) distributions should be specified for
            % specific kcat sources
            obj.params.bayesian.kcatSources         = {'dlkcat','brenda','custom'}; 
            obj.params.bayesian.sigma0logSource     = [0.4; 0.2; 0.1];
            % Shrinkage threshold: how many standard deviations a kcat should
            % deviate from its prior value to be fully updated to the posterior.
            % Higher values will make kcats "stickier" to their prior value.
            % This controls the update rate.
            obj.params.bayesian.shrinkThrDefault    = 1.5; % Default value
            obj.params.bayesian.shrinkThrSource     = [1.5, 3.5, 5.5]; % Matches order in kcatSources
            % Maximum ratio of posterior-to-prior standard deviation allowed,
            % lower values keep kcats closer to their original value. This
            % prevents random walk.
            obj.params.bayesian.varianceCapDefault  = 10;
            obj.params.bayesian.varianceCapSource   = [10,4,2];
            % If deviation < this, snap to EXACT prior. Keep certain values
            % at their starting value unless they have large impact.
            obj.params.bayesian.forcePriorThrDefault = -1; % Never force
            obj.params.bayesian.forcePriorThrSource  = [-1, 4, 8];

            % Sparsity threshold - if parameter did not move
            % much, keep it at prior. Default .5σ, alternatives
            % 0.3 (more sparse) to 0.7 (less sparse).
            obj.params.bayesian.sparsityThreshold   = 0.3;

            % Number of samples per generation
            obj.params.bayesian.scheduleGenerations = [1, 2, 9, 15];        % Schedule by which generation the sample number and target should be changed
            obj.params.bayesian.scheduleSamples     = [500, 400, 300, 200]; % Sample counts numbers corresponding to scheduleGenerations

            % Which sampled models should be selected
            obj.params.bayesian.targetAccept        = 10;  % RMSE percentile threshold (epsilon) for ABC acceptance
            obj.params.bayesian.minKeep             = 0.3; % Min fraction of samples kept each generation
            obj.params.bayesian.maxKeep             = 0.6; % Max fraction of samples kept each generation

            % Halting criteria
            obj.params.bayesian.rmseThreshold       = 0.2; % Stop when RMSE reaches this level% RMSE threshold to halt and output best posterior kcats
            obj.params.bayesian.maxGenerations      = 150;  % Hard cap on the number of ABC–SMC generations% Maximum number of generations before returning best posterior kcats
   end

        function ecModel = makeModelAnaerobic(obj,ecModel)
            % Taken from yeast-GEM 9.0.2
            ecModel = anaerobicModel_GECKO(ecModel);
        end
	    function ecModel = changeProteinBiomass(obj,ecModel,Ptot)
            % Taken from yeast-GEM 9.0.2
            ecModel = scaleBioMass_GECKO(ecModel,'protein',Ptot,'carbohydrate',false);
        end
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end
