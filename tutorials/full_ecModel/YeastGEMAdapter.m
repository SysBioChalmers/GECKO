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
            % Define initial kcat distributions (kcat * initSDmultiplDef = SD)
            obj.params.bayesian.sigma0logDefault = 2;                  % Default initial stdev of the kcat log-normal distribution
            obj.params.bayesian.kcatSources         = {'brenda','dlkcat'};  % List of annotation sources with custom SD multipliers
            obj.params.bayesian.sigma0logSelect        = [0.5; 1];          % Alterantive initial stdev of the kcat log-normal distribution

            % Number of samples per generation
            obj.params.bayesian.scheduleGenerations = [1, 2, 9, 15];         % Schedule by which generation the sample number and target should be changed
            obj.params.bayesian.scheduleSamples     = [2, 2, 100, 100]; % Sample counts numbers corresponding to scheduleGenerations

            % Which sampled models should be selected
            obj.params.bayesian.targetAccept        = 10;  % RMSE percentile threshold (epsilon) for ABC acceptance
            obj.params.bayesian.minKeep             = 0.3; % Min fraction of samples kept each generation
            obj.params.bayesian.maxKeep             = 0.6; % Max fraction of samples kept each generation

            % Low-rank proposal sampling parameters
            obj.params.bayesian.alpha               = 0.7; % Mixture weight: exploit vs explore proposal
            obj.params.bayesian.cExpl               = 3;   % Exploration inflation factor
            obj.params.bayesian.freezeStage         = 4;   % After this gen, proposal magnitudes (stds) stop adapting
            obj.params.bayesian.sigmaFloorFrac      = 0.1; % Smallest allowed std fraction relative to initial
            obj.params.bayesian.adaptFracEarly      = 0.5; % Blend factor for early adaptation of marginal scales
            obj.params.bayesian.rMax                = 150; % Maximum PCA rank in low‑rank proposal
            obj.params.bayesian.tauResidual         = 0.1; % Residual isotropic noise outside low‑rank space
            
            % Halting criteria
            obj.params.bayesian.rmseThreshold       = 0.2; % Stop when RMSE reaches this level% RMSE threshold to halt and output best posterior kcats
            obj.params.bayesian.maxGenerations      = 50;  % Hard cap on the number of ABC–SMC generations% Maximum number of generations before returning best posterior kcats
   end

        function ecModel = makeModelAnaerobic(obj,ecModel)
            % Taken from yeast-GEM 9.0.2
            ecModel = anaerobicModel_GECKO(ecModel);
        end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end
