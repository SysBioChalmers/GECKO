function [model, PdiffEnz, f]  = updateProtPool(model, measuredProt, Ptot, modelAdapter)
% updateProtPool
%   Updates the protein pool to compensate for measured proteomics data (in
%   model.ec.concs).
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   measuredProt    average sum of total measured protein in proteomics.tsv
%                   in mg/gDW, reported as protData.measuredProt by
%                   loadProtData
%   Ptot            total protein content, overwrites modelAdapter value.
%                   If specified, condition-specific fluxData.Ptot from
%                   loadFluxData can be used.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
% Output:
%   model           an ecModel where model.ec.concs is populated with
%                   protein concentrations
%   PdiffEnz        non-measured enzyme content, in mg/gDCW
%   f               f-factor as determined from proteomics data
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
         
% Ptot      the "real" total protein content in the cell
% Pmeas     the measured protein content (sum of proteomics data). Is lower
%           than Ptot, as not all protein are measured in a proteomics
%           experiment. Equals protData.measuredProt from loadProtData.
% f         the ratio enzyme/protein: how many of proteins are enzymes
% m         the ratio measured/total protein: how much is measured
% PtotEnz   the "real" total enzyme content in the cell. Equals f * Ptot
% PmeasEnz  the measured enzyme content. Equals f * Pmeas, also equals
%           sum(ecModel.ec.concs)
% PdiffEnz  the non-measured enzyme content. Equals PtotEnz - PmeasEnz.
%
% Without proteomics, protein pool is constraint by PtotEnz (= f * Ptot).
% With proteomics, protein pool should be constraint by PdiffEnz.
% [And the above two constraints are both multiplied by sigma-factor, which
% is unrelated to the adjustments made here]

% Convert ptot to mg/gDW
Ptot = Ptot * 1000;

% Calculate the protein in the model
PmeasEnz = sum(model.ec.concs,'omitnan');

if PmeasEnz == 0
    error('The model have not been constrained with proteomics')
end

% New f factor
f = PmeasEnz / measuredProt;
% Total enzyme when extrapolating from this dataset
PtotEnz = Ptot * f;

% Draw the protein in model from the measured proteins
PdiffEnz = PtotEnz - PmeasEnz;

if PdiffEnz > 0
    PdiffEnz = PdiffEnz / 1000; % Convert to g enzyme/gDCW
    Pdiff = PdiffEnz / f; % Convert to g/gDCW
    model = setProtPoolSize(model, Pdiff, f, params.sigma, modelAdapter);
    sol = solveLP(model);
    if isempty(sol.x)
        % Try relaxing sigma-factor
        model = setProtPoolSize(model, Pdiff, f, 1, modelAdapter);
        if isempty(sol.x)
            error(['Changing protein pool to ' num2str(Pdiff * f * params.sigma) ' results in a non-functional model'])
        else
            fprintf(['Relaxing of sigma-factor was required to yield a functional model.\n' ...
                     'The sigma-factor is now set to 1, run ''sigmaFitter'' using the PdiffEnz\n' ...
                     'and f outputs from updateProtPool to reduce the sigma-factor.'])
        end
    end
else
    error('The total measured protein mass exceeds the total protein content.')
end
PdiffEnz = PdiffEnz / 1000;
end
