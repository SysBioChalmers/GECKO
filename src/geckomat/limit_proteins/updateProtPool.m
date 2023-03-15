function [model, newPtot]  = updateProtPool(model, measuredProt, Ptot, modelAdapter)
% updateProtPool
%   Updates the protein pool to compensate for measured proteomics data (in
%   model.ec.concs).
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   measuredProt    average sum of total protein in proteomics.tsv in mg/gDW
%   Ptot            total protein content, overwrites modelAdapter value
%   modelAdapter    a loaded model adapter (Optional,
%
% Output:
%   model           an ecModel where model.ec.concs is populated with
%                   protein concentrations
%   newPtot         new Ptot
% Usage:
%   [model, newPtot]  = updateProtPool(model, measuredProt, Ptot, modelAdapter)

if nargin < 4 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin < 3 || isempty(Ptot)
    Ptot = params.Ptot;
    disp(['Total protein content used: ' num2str(Ptot) ' [g protein/gDw]'])
end

if nargin < 2 || isempty(measuredProt)
    % assumes not measured protein
    measuredProt = Ptot * params.f * 1000;
    disp(['Measured protein assumed: ' num2str(measuredProt/1000) ' [g protein/gDw]'])
end

% Convert ptot to mg/gDW
PtotEnz = Ptot * 1000;

% New f factor
f = measuredProt / PtotEnz;

% Calculate the protein in the model
PmeasEnz = sum(model.ec.concs,'omitnan');

if PmeasEnz == 0
    error('The model have not been constrained with proteomics')
end

% Draw the protein in model from the measured proteins
PdiffEnz = measuredProt - PmeasEnz;

if PdiffEnz > 0
    Pdiff = PdiffEnz / 1000; % Convert back to g protein/gDCW
    model = setProtPoolSize(model, Pdiff, f, params.sigma, modelAdapter);
    sol = solveLP(model);
    if isempty(sol.x)
        error(['Changing protein pool to ' num2str(Pdiff * f, params.sigma) ' resuls in a non-functional model'])
    end
else
    error('The total measured protein mass exceeds the total protein content.')
end
newPtot = PdiffEnz / 1000;
end
