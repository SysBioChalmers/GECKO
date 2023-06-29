classdef TestGEMAdapter < ModelAdapter 
    methods
        function obj = TestGEMAdapter()
            obj.params.path = fullfile(findGECKOroot,'test','unit_tests','ecTestGEM');

			obj.params.convGEM = fullfile(obj.params.path,'models','testModel.xml');

			obj.params.sigma = 0.5;

			obj.params.Ptot = 0.5;

			obj.params.gR_exp = 0.41;

			obj.params.org_name = 'testus testus';
            
            obj.params.complex.org_name = 'Testus testus';

			obj.params.kegg.ID = 'tst'; %This will not work and will not be used

            obj.params.kegg.geneID = 'kegg';

            obj.params.uniprot.type = 'taxonomy'; % 'proteome' or 'taxonomy' - will not be used
            
			obj.params.uniprot.ID = 'TE00000000'; %This will not work and will not be used
            
            obj.params.uniprot.geneIDfield = 'gene_tst';%This will not work and will not be used

            obj.params.uniprot.reviewed = true;			

			obj.params.c_source = 'E1'; 

			obj.params.bioRxn = 'R4'; %Not relevant

			obj.params.enzyme_comp = 'c';

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
