%Abstract Base class for adapters for different species
classdef (Abstract) ModelAdapter 
	methods (Abstract)
        result = getFilePath(obj, filename)
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
        
		function model = manualModifications(obj,model) %default is to do nothing
        end
        
        %The genes returned here should match the gene id set in the parameter params.uniprotGeneIdField
        function genes = getUniprotCompatibleGenes(obj,inGenes)
            genes = inGenes;
        end
        
        function folder = getBrendaDBFolder(obj)
            folder = fullfile(findGECKOroot(),'databases');
        end
        
        function x = getPhylDistStructPath(obj)
            x =  fullfile(findGECKOroot(),'databases','PhylDist.mat');
        end
        
    end

    %To have the params public is a bit "ugly", but very practical 
    %if we want to change a parameter
    properties (Access = public) 
        params;
    end
end
