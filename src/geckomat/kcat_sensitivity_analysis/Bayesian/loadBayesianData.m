function bayData = loadBayesianData(varargin)
% loadBayesianData  Load experimental data for Bayesian kcat tuning.
%
% Reads the experimental data files used by bayesianSensitivityTuning from
% the data subfolder of the model directory defined by the model adapter.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% bayData : struct
%     structure with the loaded experimental data.
%
% Notes
% -----
% The bayData structure has the following fields:
%
% - fluxData : flux data loaded from bayesianFluxData.tsv, with the biomass
%   reaction id stored in fluxData.biomass.
% - maxGrate : maximum growth-rate data loaded from bayesianMaxGrowth.tsv,
%   with the biomass reaction id stored in maxGrate.biomass.
% - zeroFlux : exchange reactions assumed to carry zero flux, loaded from
%   bayesianZeroExch.tsv.
%
% See also
% --------
% bayesianSensitivityTuning, abc_max

p = parseGECKOargs(varargin, { ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

basePath            = modelAdapter.params.path;
bayData.fluxData = []; bayData.maxGrate = []; bayData.zeroFlux = [];
bayData.fluxData = loadFluxData(fullfile(basePath,'data','bayesianFluxData.tsv'));
if ~isempty(bayData.fluxData)
    bayData.fluxData.biomass = modelAdapter.params.bioRxn;
end
bayData.maxGrate = loadFluxData(fullfile(basePath,'data','bayesianMaxGrowth.tsv'));
bayData.maxGrate.biomass = modelAdapter.params.bioRxn;
bayData.zeroFlux = table2cell(readtable(fullfile(basePath,'data','bayesianZeroExch.tsv'), 'Delimiter', '\t', 'FileType','delimitedtext'));
end
