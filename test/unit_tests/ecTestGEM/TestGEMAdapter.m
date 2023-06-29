classdef TestGEMAdapter < ModelAdapter 
    methods
        function obj = TestGEMAdapter()
            %Set initial values of the obj.params - they can be changed by the user
            
            %Directory where all model-specific files and scripts are kept.
            %Is assumed to follow the GECKO-defined folder structure. The
            %code below refers to tutorials/tutorial_yeast-GEM in the GECKO path.
            obj.params.path = fullfile(findGECKOroot,'test','unit_tests','ecTestGEM');

			%Path to the conventional GEM that this ecModel will be based on.
			obj.params.convGEM = fullfile(obj.params.path,'models','testModel.xml');

            %Average enzyme saturation factor
			obj.params.sigma = 0.5;

			%Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.5;      %Assumed constant

			%Minimum growth rate the model should grow at [1/h]
			obj.params.gR_exp = 0.41;     %[g/gDw h] 

			%Provide your organism scientific name
			obj.params.org_name = 'testus testus';
            
            %Matching name for Complex Portal
            obj.params.complex.org_name = 'Testus testus';

			%Provide your organism KEGG ID, selected at
			%https://www.genome.jp/kegg/catalog/org_list.html
			obj.params.kegg.ID = 'tst'; %This will not work and will not be used
            obj.params.kegg.geneID = 'kegg';

			%Provide what identifier should be used to query UniProt.
            %Select proteome IDs at https://www.uniprot.org/proteomes/
            %or taxonomy IDs at https://www.uniprot.org/taxonomy.
            obj.params.uniprot.type = 'taxonomy_id'; % 'proteome' or 'taxonomy_id' - will not be used
			obj.params.uniprot.ID = 'TE00000000'; %This will not work and will not be used
            
            %Field for Uniprot gene id - should match the gene ids used in the 
            %GPRs. Note that this is a field in the web request to uniprot - 
            %it has to match one of the fields there
            obj.params.uniprot.geneIDfield = 'gene_tst';%This will not work and will not be used

            %Whether only reviewed data from UniProt should be considered.
            %Reviewed data has highest confidence, but coverage might be (very)
            %low for non-model organisms
            obj.params.uniprot.reviewed = true;			

			%The name of the exchange reaction that supplies the model with carbon (rxnNames)
			obj.params.c_source = 'E1'; 

			%Rxn Id for biomass pseudoreaction
			obj.params.bioRxn = 'R4'; %Not relevant

			%Compartment name in which the added enzymes should be located
			obj.params.enzyme_comp = 'c';

            %The pool size, fitted to data
            obj.params.f = 4;

        end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = false(length(model.rxns), 1);
			spontRxnNames = rxns_tsv.rxns;
			spont(5) = true; %R4 is spontaneous
        end
        
        function folder = getBrendaDBFolder(obj)
            folder = fullfile(obj.params.path,'data');
        end
        
        function x = getPhylDistStructPath(obj)
            x =  fullfile(obj.params.path,'data','PhylDist.mat');
        end

	end
end
