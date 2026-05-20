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

            %% OpenKineticsPredictor settings (used by submit/fetchOpenKineticsPredictor)
            % These are optional; the functions fall back to these same
            % defaults if the okp block is absent. The API key is NOT set
            % here (it is a secret): provide it as a function argument, the
            % OKP_API_KEY environment variable, or data/okpApiKey.txt.
            % Predictor method. One of: CataPro, CatPred, DLKcat, EITLEM,
            % KinForm-H, KinForm-L, UniKP (see GET /api/v1/methods/).
            obj.params.okp.method                   = 'CataPro';
            % Kinetic parameter(s) to predict. geckopy/GECKO use kcat.
            obj.params.okp.targets                  = {'kcat'};
            % How to handle sequences exceeding a method's max length.
            obj.params.okp.handleLongSequences      = 'truncate';
            % Append per-row similarity-to-training-data columns to output.
            obj.params.okp.includeSimilarityColumns = true;
            % Canonicalize substrate SMILES server-side before prediction.
            obj.params.okp.canonicalizeSubstrates   = true;

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