% KEY_CLASSNAME  Template ModelAdapter to copy and adapt per organism.
%
% Starting point for a species-specific adapter. Copy this file, rename
% the class, and edit the parameter values set in the constructor and the
% getSpontaneousReactions method to match the target organism and model.
%
% See also
% --------
% ModelAdapter, ModelAdapterManager
classdef KEY_CLASSNAME < ModelAdapter
    methods
        function obj = KEY_CLASSNAME()
            % KEY_CLASSNAME  Construct the adapter and set default parameters.
            %
            % Sets initial values of obj.params; they can be changed by
            % the user after construction.
            %
            % Returns
            % -------
            % obj : ModelAdapter
            %     the constructed adapter instance.

            % Directory where all model-specific files and scripts are kept.
            % Is assumed to follow the GECKO-defined folder structure.
            obj.params.path = fullfile('KEY_PATH', 'KEY_NAME');
            addpath(fullfile(obj.params.path,'code'));

			% Path to the conventional GEM that this ecModel will be based on.
			obj.params.convGEM = fullfile(obj.params.path,'models','yourModel.xml');

			% Average enzyme saturation factor
			obj.params.sigma = 0.5;

			% Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.5;

			% Fraction of enzymes in the model [g enzyme/g protein]
			obj.params.f = 0.5;
            
            % Growth rate the model should be able to reach when not
            % constraint by nutrient uptake (e.g. max growth rate) [1/h]
			obj.params.gR_exp = 0.41;

			% Provide your organism scientific name
			obj.params.org_name = 'genus species';
            
            % Taxonomic identifier for Complex Portal
            obj.params.complex.taxonomicID = [];

			% Provide your organism KEGG ID, selected at
			% https://www.genome.jp/kegg/catalog/org_list.html
			obj.params.kegg.ID = 'sce';
            % Field for KEGG gene identifier; should match the gene
            % identifiers used in the model. With 'kegg', it takes the
            % default KEGG Entry identifier (for example YER023W here:
            % https://www.genome.jp/dbget-bin/www_bget?sce:YER023W).
            % Alternatively, gene identifiers from the "Other DBs" section
            % of the KEGG page can be selected. For example "NCBI-GeneID",
            % "UniProt", or "Ensembl". Not all DB entries are available for
            % all organisms and/or genes.
            obj.params.kegg.geneID = 'kegg';

			% Provide what identifier should be used to query UniProt.
            % Select proteome IDs at https://www.uniprot.org/proteomes/
            % or taxonomy IDs at https://www.uniprot.org/taxonomy.
            obj.params.uniprot.type = 'taxonomy'; % 'proteome' or 'taxonomy'
			obj.params.uniprot.ID = '559292'; % should match the ID type
            % Field for Uniprot gene ID - should match the gene ids used in the 
            % model. It should be one of the "Returned Field" entries under
            % "Names & Taxonomy" at this page: https://www.uniprot.org/help/return_fields
            obj.params.uniprot.geneIDfield = 'gene_oln';
            % Whether only reviewed data from UniProt should be considered.
            % Reviewed data has highest confidence, but coverage might be (very)
            % low for non-model organisms
            obj.params.uniprot.reviewed = false;

			% Reaction ID for glucose exchange reaction (or other preferred carbon source)
			obj.params.c_source = 'r_1714'; 

			% Reaction ID for biomass pseudoreaction
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
            % (rmseFlux + w*rmseMaxGrowth) / (w+1). At 1 both count equally.
            obj.params.evotune.maxGrowthWeight     = 1.0;
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
            obj.params.evotune.maxGenerations      = 150;
        end

        % function ecModel = makeModelAnaerobic(ecModel)
        %     % Define a model-specific function in the 'code' subfolder,
        %     % that can constrain the model to anaerobic conditions, and
        %     % include the name of this function here. This is used by
        %     % cmaesKcatTuning.m (via evotuneScore.m) if the fluxData
        %     % has anaerobic conditions.
        %     addpath(fullfile(obj.params.path,'code'));
        %     ecModel = nameOfModelSpecificFunction(ecModel);
        % end

        function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
            % getSpontaneousReactions  Identify spontaneous reactions in the model.
            %
            % Indicates how spontaneous reactions are identified. Here it
            % is done by the reaction having 'spontaneous' in its name.
            %
            % Parameters
            % ----------
            % model : struct
            %     a model in RAVEN format.
            %
            % Returns
            % -------
            % spont : logical
            %     true for each reaction identified as spontaneous.
            % spontRxnNames : cell
            %     names of the reactions identified as spontaneous.
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end