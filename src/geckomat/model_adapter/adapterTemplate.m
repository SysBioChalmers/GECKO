classdef KEY_CLASSNAME < ModelAdapter
    methods
        function obj = KEY_CLASSNAME()
            % Set initial values of the obj.params - they can be changed by the user
            
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

            %% Hyperparameters for Bayesian kcat fitting
            % Define initial kcat distributions (kcat * initSDmultiplDef = SD)
            obj.params.bayesian.initSDmultiplDef    = 1;                    % Default initial SD 
            obj.params.bayesian.kcatSources         = {'brenda','dlkcat'};  % List of annotation sources with custom SD multipliers
            obj.params.bayesian.initSDmultipl       = [0.05; 0.3];          % Multipliers that overwrite initSDmultplDef, matching kcatSources

            % Number of samples per generation
            obj.params.bayesian.scheduleGenerations = [1, 2, 9, 15];         % Schedule by which generation the sample number and target should be changed
            obj.params.bayesian.scheduleSamples     = [1500, 500, 300, 200]; % Sample counts numbers corresponding to scheduleGenerations

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

        % function ecModel = makeModelAnaerobic(ecModel)
        %     % Define a model-specific function in the 'code' subfolder,
        %     % that can constrain the model to anaerobic conditions, and
        %     % include the name of this function here. This is used by
        %     % bayesianSensitivityTuning.m (via abc_max.m) if the fluxData
        %     % has anaerobic conditions.
        %     addpath(fullfile(obj.params.path,'code'));
        %     ecModel = nameOfModelSpecificFunction(ecModel);
        % end

        function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
            % Indicates how spontaneous reactions are identified. Here it
            % is done by the reaction have 'spontaneous' in its name.
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end