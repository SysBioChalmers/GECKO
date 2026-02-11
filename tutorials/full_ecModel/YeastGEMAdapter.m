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
            obj.params.bayesian.scheduleGenerations = [1, 2, 9, 15]; % Schedule by which generation the sample number and target should be changed
            obj.params.bayesian.scheduleSamples     = [1500, 500, 300, 200]; % Schedule of sample numbers (matching scheduleGenerations)
            obj.params.bayesian.scheduleTarget      = [0.4, 0.2, 0.1, 0.05]; % Schedule of target acceptance fractions (matching scheduleGenerations)
            obj.params.bayesian.minKeep             = 10; % Minimum number of samples to keep
            obj.params.bayesian.rmseThreshold       = 0.2; % RMSE threshold to halt and output best posterior kcats
            obj.params.bayesian.maxGenerations      = 200; % Maximum number of generations before returning best posterior kcats
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
