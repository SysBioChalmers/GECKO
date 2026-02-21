function bayData = loadBayesianData(modelAdapter)

if nargin < 1 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

basePath            = modelAdapter.params.path;
bayData.fluxData = loadFluxData(fullfile(basePath,'data','bayesianFluxData.tsv'));
bayData.fluxData.biomass = modelAdapter.params.bioRxn;
bayData.maxGrate = loadFluxData(fullfile(basePath,'data','bayesianMaxGrowth.tsv'));
bayData.maxGrate.biomass = modelAdapter.params.bioRxn;
bayData.zeroFlux = table2cell(readtable(fullfile(basePath,'data','bayesianZeroExch.tsv'), 'Delimiter', '\t', 'FileType','delimitedtext'));
end
