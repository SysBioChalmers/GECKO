function [fluxData, maxGrate, zeroFlux] = loadBayesianData(modelAdapter)

if nargin < 1 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

basePath            = modelAdapter.params.path;
fluxData = loadFluxData(fullfile(basePath,'data','bayesianFluxData.tsv'));
fluxData.biomass = modelAdapter.params.bioRxn;
maxGrate = loadFluxData(fullfile(basePath,'data','bayesianMaxGrowth.tsv'));
maxGrate.biomass = modelAdapter.params.bioRxn;
zeroFlux = table2cell(readtable(fullfile(basePath,'data','bayesianZeroExch.tsv'), 'Delimiter', '\t', 'FileType','delimitedtext'));
end
