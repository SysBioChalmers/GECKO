function parameters = getModelParameters
% getModelParameters
%
%   Set model and organism specific parameters that are used by the
%   ecModel generation pipeline.
%
%   Ivan Domenzain. Last edited: 2019-07-12
%

%Average enzyme saturation factor
parameters.sigma    = 0.5;
%Total protein content in the cell [g protein/gDw]
parameters.Ptot     = 0.5;      %Assumed constant
%Minimum growth rate the model should grow at [1/h]
parameters.gR_exp   = 0.41;     %[g/gDw h] 
%Provide your organism scientific name
parameters.org_name = 'saccharomyces cerevisiae';
%Provide your organism KEGG ID
parameters.keggID   = 'sce';
%The name of the exchange reaction that supplies the model with carbon
%(rxnNames)
parameters.c_source = 'D-glucose exchange (reversible)'; 
%Rxn name for biomass pseudoreaction
parameters.bioRxn   = '';
%Compartment name in which the added enzymes should be located
parameters.enzyme_comp = '';

end