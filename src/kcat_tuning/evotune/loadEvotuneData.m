function evoData = loadEvotuneData(varargin)
% loadEvotuneData  Load experimental data for evotune (CMA-ES) kcat tuning.
%
% Reads the three experimental data files used by cmaesKcatTuning from the
% data subfolder of the model directory defined by the model adapter. Each
% file is optional: a missing file yields [] (fluxData/maxGrate) or an
% empty cell array (zeroFlux) rather than an error.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% evoData : struct
%     structure with the loaded experimental data.
%
% Notes
% -----
% The evoData structure has the following fields:
%
% - fluxData : flux data loaded from evotuneFluxData.tsv, with the biomass
%   reaction id stored in fluxData.biomass. [] if that file isn't present.
% - maxGrate : maximum growth-rate data loaded from evotuneMaxGrowth.tsv,
%   with the biomass reaction id stored in maxGrate.biomass. [] if that
%   file isn't present.
% - zeroFlux : reaction IDs assumed to carry zero flux, loaded from
%   evotuneZeroExch.tsv. Empty cell array if that file isn't present.
%
% See also
% --------
% cmaesKcatTuning, evotuneScore
p = parseGECKOargs(varargin, { ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

basePath = modelAdapter.params.path;

evoData.fluxData = loadFluxData(fullfile(basePath,'data','evotuneFluxData.tsv'), modelAdapter);
if ~isempty(evoData.fluxData)
    evoData.fluxData.biomass = modelAdapter.params.bioRxn;
end

evoData.maxGrate = loadFluxData(fullfile(basePath,'data','evotuneMaxGrowth.tsv'), modelAdapter);
if ~isempty(evoData.maxGrate)
    evoData.maxGrate.biomass = modelAdapter.params.bioRxn;
end

% Unlike loadFluxData (used for fluxData/maxGrate above, which already
% degrades to [] on a missing file), readtable has no such guard, so this
% is wrapped explicitly to keep all three files uniformly optional.
zeroExchFile = fullfile(basePath,'data','evotuneZeroExch.tsv');
if isfile(zeroExchFile)
    evoData.zeroFlux = table2cell(readtable(zeroExchFile, 'Delimiter', '\t', 'FileType','delimitedtext'));
else
    evoData.zeroFlux = {};
end
end
