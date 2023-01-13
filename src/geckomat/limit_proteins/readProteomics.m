function model = readProteomics(model,protData)
% readProteomics
%   Reads absolute proteomics data, from either a file or protData
%   structure. The proteins should be annotated with Uniprot IDs, while the
%   protein levels should be given as mmol/gDCW. Any existing protein
%   concentrations in model.ec.concs are discarded. If no data is provided
%   a particular protein, its level is NaN.
%
% Input:
%   model       an ec-model
%   protData    either string where to find the input file, or a structure
%               with the following fields:
%               protData.uniprot  Uniprot identifiers
%               protData.level    Protein levels in mmol/gDCW
%
% Output:
%   model       an ec-model where model.ec.concs is populated with protein
%               concentrations.
%
% Note: to also constrain the model with the content of model.ec.concs, you
% should run constrainProtConcs.

% Here some implementation of ModelAdapter, to find default userData
% folder.

switch class(protData)
    case {'char','string'}
        fID = fopen(fullfile(protData),'r');
        fileContent = textscan(fID, '%s %s', 'HeaderLines', 1);
        fclose(fID);
        uniprot = fileContent{1};
        level   = str2double(fileContent{2});
    case 'struct'
        uniprot = protData.uniprot;
        level   = protData.level;
    otherwise
        error(['protData should either be a structure or string, not a ' class(protData)])
end

%Redefine an empty model.ec.concs vector
model.ec.concs=nan(numel(model.ec.enzymes),1);

[a,b] = ismember(uniprot, model.ec.enzymes);
model.ec.concs(b(a)) = level(a);
end