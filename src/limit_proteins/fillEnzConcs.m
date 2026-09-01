function model = fillEnzConcs(model, protData, varargin)
% fillEnzConcs  Populate model.ec.concs from proteome data.
%
% Uses the protein concentrations from protData to fill model.ec.concs.
% Protein levels should be provided in mg/gDCW. If no data is provided for a
% particular protein, its level is NaN. Existing entries in model.ec.concs
% are overwritten.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% protData : struct
%     structure with proteome data, from loadProtFluxData, with the fields:
%
%     - uniprotIDs : cell array with Uniprot IDs matching protData.abundances.
%     - abundances : matrix of proteomics data.
%
% Name-Value Arguments
% --------------------
% dataCol : double
%     number indicating the column in protData.abundances that contains the
%     relevant protein concentrations. protData may contain data from
%     multiple conditions/samples/experiments, each with their own column in
%     protData.abundances (default 1).
%
% Returns
% -------
% model : struct
%     an ecModel where model.ec.concs is populated with protein
%     concentrations.
%
% Notes
% -----
% To also constrain the model with the content of model.ec.concs, you should
% run constrainEnzConcs.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = fillEnzConcs(model, protData);
%     model = fillEnzConcs(model, protData, 'dataCol', 2);
%
% See also
% --------
% constrainEnzConcs, loadProtFluxData

p = parseGECKOargs(varargin, { ...
    'dataCol', 1});
dataCol = p.dataCol;
if isempty(dataCol); dataCol = 1; end

uniprotIDs = protData.uniprotIDs;
abundances = protData.abundances(:,dataCol);

%Redefine an empty model.ec.concs vector
model.ec.concs=nan(numel(model.ec.enzymes),1);

[a,b] = ismember(uniprotIDs, model.ec.enzymes);
model.ec.concs(b(a)) = abundances(a);
end
