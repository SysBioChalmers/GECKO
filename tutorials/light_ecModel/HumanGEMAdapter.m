classdef HumanGEMAdapter < ModelAdapter 
    methods
        function obj = HumanGEMAdapter()
            geckoPath = findGECKOroot;
            obj.params.path = fullfile(geckoPath,'tutorials','light_ecModel');

            % The model distributed with the light_ecModel tutorial is Human-GEM
            % is version 1.15.0, available from
            % https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.3.0
            % In addition, the following lines were run to reduce its size
            % before storing it in this GECKO tutorial:
            % ihuman = simplifyModel(ihuman,false,false,true,true);
            % [ihuman.grRules,skipped] = simplifyGrRules(ihuman.grRules,true);
            % ihuman = deleteUnusedGenes(ihuman);
            % ihuman = rmfield(ihuman,{'unconstrained','rxnReferences','rxnFrom','metFrom','rxnConfidenceScores'});
            % ihuman.name = ihuman.id;
            % writeYAMLmodel(ihuman,'human-GEM.yml')
            obj.params.convGEM = fullfile(obj.params.path,'models','human-GEM.yml');

			obj.params.sigma = 0.1; %This was changed to a low number to give a reasonable growth rate - this should be investigated more

			obj.params.Ptot = 0.505717472;  %Average across NCI60 cell lines

			obj.params.f = 0.412; %estimated as TPM of model genes / all genes in GTEx

			obj.params.gR_exp = 0.020663429; %[g/gDw h]/Average across NCI60 cell lines
            
			obj.params.c_UptakeExp = 0.641339301; %[mmol/gDw h]/Average across NCI60 cell lines
            
			obj.params.org_name = 'homo sapiens';
            
            obj.params.complex_org_name = 'Homo sapiens';

			obj.params.keggID = 'hsa';

            obj.params.keggGeneIdentifier = 'Ensembl';

            obj.params.uniprotIDtype = 'proteome';

			obj.params.uniprotID = 'UP000005640';

            obj.params.uniprotGeneIdField = 'gene_primary';

            obj.params.uniprotReviewed = true;
            
			obj.params.c_source = 'MAR09034'; 

			obj.params.bioRxn = 'MAR13082';

			obj.params.NGAM = 'MAR09931'; %This is not used

			obj.params.enzyme_comp = 'Cytosol';
            
        end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
            fID=fopen(fullfile(obj.params.path,'data','spontaneousReactions.tsv'),'r');
            rxns_tsv = textscan(fID,'%q %q','Delimiter','\t','HeaderLines',1);
            fclose(fID);
            % The above file is derived from the reactions.tsv that is available
            % in the Human-GEM repository, matching version 1.15.0 of the model.
            % https://github.com/SysBioChalmers/Human-GEM/blob/194ebe5431c83e25f78df61caacad2fa485b5cb4/model/reactions.tsv
            % Only the columns with reaction identifiers and assignment of
            % spontaneous reactions were kept in spontaneousReactions.tsv.
            spont = logical(str2double(cell2mat(rxns_tsv{2})));
			spontRxnNames = rxns_tsv{1};
            [spont,~] = ismember(model.rxns,spontRxnNames(spont));
            spontRxnNames = model.rxns(spont);
        end
    end
    
    methods(Static)
        function path = getHumanGEMRootPath()
            path = fileparts(which('Human-GEM.mat'));% This file should be on the path
            path = strrep(path, '\', '/'); % get rid of backslashes in Windows
            if ~endsWith(path, '/')
                path = strcat(path,'/');
            end
            % Now remove the model/ at the end
            path = (path(1:strlength(path)-6));
        end
	end
end
