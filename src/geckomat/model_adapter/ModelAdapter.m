%Abstract Base class for adapters for different species
classdef (Abstract) ModelAdapter
    methods (Abstract)
        [spont,spontRxnNames] = getSpontaneousReactions(obj,model);
    end
    methods
        %not really needed, we could access it directly. Not very nice, but practical
        function parameters = getParameters(obj)
            parameters = obj.params;
        end

        %This allows for changing parameters in an adapter without
        %creating a new class. May be convenient, although less clean
        %function obj = setParameters(obj, parameters)
        %    obj.mParams = parameters;
        %end

        %The genes returned here should match the gene id set in the parameter params.uniprot.geneIDfield
        function genes = getUniprotCompatibleGenes(obj,inGenes)
            genes = inGenes;
        end

        function uniprotIDs = getUniprotIDsFromTable(obj,modelGenes)
            conversionTable = fullfile(obj.params.path,'data','uniprotConversion.tsv');
            if exist(conversionTable,'file')
                fID=fopen(conversionTable,'r');
                conversionTable = textscan(fID,'%q %q','Delimiter','\t','HeaderLines',1);
                fclose(fID);

                modelIDs    = conversionTable{1,1};
                uniprots    = conversionTable{1,2};

                [a,b] = ismember(modelGenes,modelIDs);

                uniprotIDs  = strings(numel(modelGenes),1);
                uniprotIDs(a)  = uniprots(b(a));
                disp('The model genes are matched to Uniprot via the table at data/uniprotConversion.tsv.')
            else
                uniprotIDs = modelGenes;
            end
        end

        function folder = getBrendaDBFolder(obj)
            folder = fullfile(findGECKOroot(),'databases');
        end

        function x = getPhylDistStructPath(obj)
            x =  fullfile(findRAVENroot(),'external','kegg','keggPhylDist.mat');
        end
    end

    %To have the params public is a bit "ugly", but very practical
    %if we want to change a parameter
    properties (Access = public)
        params;
    end
end
