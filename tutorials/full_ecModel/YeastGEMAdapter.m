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

            %% Hyperparameters for evotune (CMA-ES) kcat fitting
            % Default initial uncertainty (standard deviation in log-space) for kcat values
            obj.params.evotune.sigma0logDefault    = 0.5;
            % Data sources for kcat values, ordered from least to most trusted
            obj.params.evotune.kcatSources         = {'dlkcat','brenda','custom'};
            % Initial uncertainty for each source (lower = more trusted data)
            obj.params.evotune.sigma0logSource     = [0.4; 0.2; 0.1];

            % Weight on the max-growth RMSE against the flux RMSE:
            % (rmseFlux + w*rmseMaxGrowth) / (w+1). At 2 the max-growth
            % conditions count double against the flux conditions.
            obj.params.evotune.maxGrowthWeight     = 2.0;
            % Weight on a prior term in the search objective,
            % rmse + w*mean((log(k/k0)/sigma0log)^2). Keeps large corrections
            % reproducible across seeds; 0 scores on RMSE alone.
            obj.params.evotune.priorPenaltyWeight  = 0.03;
            % Give isozyme copies of one reaction (sharing a prior value and
            % a source) a single shared kcat, so the search cannot invent a
            % distinction the kcat assignment never made.
            obj.params.evotune.tieIsozymes         = true;

            % Stop optimization when RMSE falls below this threshold
            obj.params.evotune.rmseThreshold       = 0.2;
            % Maximum number of CMA-ES generations before termination
            obj.params.evotune.maxGenerations      = 212;
   end

        function ecModel = makeModelAnaerobic(obj,ecModel)
            % Constrains ecModel to anaerobic growth conditions (blocked
            % oxygen uptake, enabled sterol/fatty acid exchanges, and other
            % curations), as implemented in anaerobicModel_GECKO.
            ecModel = anaerobicModel_GECKO(ecModel);
        end
	    function ecModel = changeProteinBiomass(obj,ecModel,Ptot)
            % Currently a no-op: returns ecModel unchanged. Intended to
            % rescale the biomass reaction's protein content to Ptot (as
            % scaleBioMass_GECKO would), but this step is disabled for now.
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
