function exportBayesianResults(sigmaLogTrace, kcatTrace, rmseTrace, diagnostics, posteriorSamples, modelAdapter, ecModel, outputPath)
% exportBayesianResults  Write bayesianSensitivityTuning outputs to flat TSV files.
%
%   Saves the traces, diagnostics and posterior ensemble produced by
%   bayesianSensitivityTuning into tab-separated text files under
%   fullfile(modelAdapter.params.path, 'output') so they can be archived or
%   loaded into any other tool (Excel, Python/pandas, R) without MATLAB.
%
% INPUTS (in this order):
%   sigmaLogTrace    [D x (nGen+1)] posterior sigmaLog per generation (col 1 = prior)
%   kcatTrace        [D x (nGen+1)] posterior kcat means per generation (col 1 = prior)
%   rmseTrace        [1 x (nGen+1)] best RMSE per generation (index 1 = prior)
%   diagnostics      struct from bayesianSensitivityTuning
%   posteriorSamples struct with fields .kcats [D x Nacc] and .rmse [Nacc x 1]
%   modelAdapter     (optional) model adapter; used for params.path. Defaults
%                    to ModelAdapterManager.getDefault().
%   ecModel          (optional) final ecModel; if given, parameter-indexed
%                    files also get reaction_id and source columns.
%   outputPath       (optional) output directory. Defaults to
%                    fullfile(modelAdapter.params.path, 'output').
%
% OUTPUT FILES (tab-separated .tsv, plus one .txt summary):
%   rmse_trace.tsv                 generation, rmse
%   kcat_trace.tsv                 param rows x generation cols (first col "prior")
%   sigmalog_trace.tsv             param rows x generation cols
%   shrinkage_trace.tsv            param rows x generation cols (diagnostics.shrinkageTrace)
%   diagnostics_pergeneration.tsv  generation x per-generation scalar metrics
%   diagnostics_bysource.tsv       long format: generation, source, per-source metrics
%   posterior_kcats.tsv            param rows x accepted-sample cols
%   posterior_rmse.tsv             sample_index, rmse
%   run_metadata.txt               key/value summary of the run
%
% NOTE on generations: the trace files index generations 0..nGen with column 1
%   being the prior. The diagnostics arrays start at generation 1 and, after an
%   early stop, are slightly shorter than the traces. Each file is written at
%   its own width and labels its generation column explicitly, so the files can
%   be joined on the integer generation number.
%
% Example:
%   [ecModel, rmseTrace, kcatTrace, sigmaLogTrace, diagnostics, posteriorSamples] = ...
%       bayesianSensitivityTuning(ecModel);
%   exportBayesianResults(sigmaLogTrace, kcatTrace, rmseTrace, diagnostics, ...
%       posteriorSamples, modelAdapter, ecModel);

%% Resolve optional arguments
if nargin < 6 || isempty(modelAdapter)
    try
        modelAdapter = ModelAdapterManager.getDefault();
    catch
        modelAdapter = [];
    end
end
if nargin < 7
    ecModel = [];
end
if nargin < 8
    outputPath = '';
end

%% Resolve output directory: fullfile(params.path,'output')
if isempty(outputPath)
    try
        outputPath = fullfile(modelAdapter.params.path, 'output');
    catch
        outputPath = fullfile(pwd, 'output');
        warning('exportBayesianResults:NoPath', ...
            'Could not resolve modelAdapter.params.path; writing to %s instead.', outputPath);
    end
end
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

%% Light input validation (the input order is unusual, so guard the obvious slip)
if ~isvector(rmseTrace)
    error('exportBayesianResults:BadInput', ...
        ['rmseTrace must be a vector. Check argument order: ', ...
         '(sigmaLogTrace, kcatTrace, rmseTrace, diagnostics, posteriorSamples, ...).']);
end
if ~isempty(kcatTrace) && ~isempty(sigmaLogTrace) && ...
        size(kcatTrace,1) ~= size(sigmaLogTrace,1)
    warning('exportBayesianResults:RowMismatch', ...
        'kcatTrace and sigmaLogTrace have different row counts (%d vs %d); writing anyway.', ...
        size(kcatTrace,1), size(sigmaLogTrace,1));
end

D = size(kcatTrace, 1);

%% Per-parameter row labels (reaction_id, source) from ecModel if available
paramIDs = {};
paramSources = {};
if ~isempty(ecModel) && isfield(ecModel, 'ec')
    if isfield(ecModel.ec, 'rxns') && numel(ecModel.ec.rxns) == D
        paramIDs = ecModel.ec.rxns;
    end
    if isfield(ecModel.ec, 'source') && numel(ecModel.ec.source) == D
        paramSources = ecModel.ec.source;
    end
end

numFmtKcat = '%.8g';   % kcat / sigma / posterior values
numFmtDiag = '%.6g';   % diagnostics (rates, counts, fractions)

filesWritten = {};

%% 1. RMSE trace  (generation 0..nGen, index 1 = prior)
rmseVec = rmseTrace(:);
genR = (0:numel(rmseVec)-1)';
if writeGenTable(fullfile(outputPath,'rmse_trace.tsv'), genR, {'rmse'}, rmseVec, numFmtKcat)
    filesWritten{end+1} = 'rmse_trace.tsv';
end

%% 2. kcat trace (params x generations, first data column = prior)
genHeadersTrace = [{'prior'}, arrayfun(@(g) sprintf('gen%d',g), ...
    1:(size(kcatTrace,2)-1), 'UniformOutput', false)];
if writeParamMatrixTSV(fullfile(outputPath,'kcat_trace.tsv'), kcatTrace, ...
        genHeadersTrace, paramIDs, paramSources, numFmtKcat)
    filesWritten{end+1} = 'kcat_trace.tsv';
end

%% 3. sigmaLog trace (params x generations)
if ~isempty(sigmaLogTrace)
    genHeadersSig = [{'prior'}, arrayfun(@(g) sprintf('gen%d',g), ...
        1:(size(sigmaLogTrace,2)-1), 'UniformOutput', false)];
    if writeParamMatrixTSV(fullfile(outputPath,'sigmalog_trace.tsv'), sigmaLogTrace, ...
            genHeadersSig, paramIDs, paramSources, numFmtKcat)
        filesWritten{end+1} = 'sigmalog_trace.tsv';
    end
end

%% 4. Per-parameter shrinkage trace (diagnostics.shrinkageTrace: params x generations)
if isfield(diagnostics, 'shrinkageTrace') && ~isempty(diagnostics.shrinkageTrace)
    nShr = size(diagnostics.shrinkageTrace, 2);
    shrHeaders = arrayfun(@(g) sprintf('gen%d',g), 1:nShr, 'UniformOutput', false);
    if writeParamMatrixTSV(fullfile(outputPath,'shrinkage_trace.tsv'), ...
            diagnostics.shrinkageTrace, shrHeaders, paramIDs, paramSources, numFmtDiag)
        filesWritten{end+1} = 'shrinkage_trace.tsv';
    end
end

%% 5. Per-generation scalar diagnostics (one row per generation)
scalarFields = {'acceptanceRateTrace','proposalAcceptRate','boundaryViolationsTrace', ...
    'directednessTrace','epsilonTrace','nSamplesTrace','nAcceptedTrace', ...
    'meanAcceptedRMSE','medianAcceptedRMSE','sparsityCountTrace', ...
    'sparsityFractionTrace','diversityTrace','proposalWidthTrace'};
nDiag = 0;
if isstruct(diagnostics)
    present = scalarFields(isfield(diagnostics, scalarFields));
    if ~isempty(present)
        nDiag = size(diagnostics.(present{1}), 2);   % common generation count
        keep = true(1, numel(present));
        dataMat = zeros(nDiag, numel(present));
        for k = 1:numel(present)
            v = diagnostics.(present{k});
            if numel(v) == nDiag
                dataMat(:,k) = v(:);
            else
                keep(k) = false;   % skip a field of unexpected length
            end
        end
        present = present(keep);
        dataMat = dataMat(:, keep);
        genD = (1:nDiag)';
        if writeGenTable(fullfile(outputPath,'diagnostics_pergeneration.tsv'), ...
                genD, present, dataMat, numFmtDiag)
            filesWritten{end+1} = 'diagnostics_pergeneration.tsv';
        end
    end
end

%% 6. Per-source diagnostics (long/tidy: generation, source, metrics)
srcFields = {'activeBySource','nearPriorBySource','meanDeviationBySource', ...
    'varianceRatioBySource','changeMagnitudeQ25','changeMagnitudeQ50','changeMagnitudeQ75'};
if isstruct(diagnostics) && isfield(diagnostics, 'sourceLabels')
    srcPresent = srcFields(isfield(diagnostics, srcFields));
    if ~isempty(srcPresent)
        sourceLabels = diagnostics.sourceLabels;
        M1 = diagnostics.(srcPresent{1});
        nG = size(M1, 2);
        nSrcRows = min(numel(sourceLabels), size(M1, 1));
        fid = fopen(fullfile(outputPath,'diagnostics_bysource.tsv'), 'w');
        if fid ~= -1
            fprintf(fid, 'generation\tsource');
            for k = 1:numel(srcPresent), fprintf(fid, '\t%s', srcPresent{k}); end
            fprintf(fid, '\n');
            for g = 1:nG
                for s = 1:nSrcRows
                    fprintf(fid, '%d\t%s', g, sourceLabels{s});
                    for k = 1:numel(srcPresent)
                        M = diagnostics.(srcPresent{k});
                        fprintf(fid, ['\t' numFmtDiag], M(s, g));
                    end
                    fprintf(fid, '\n');
                end
            end
            fclose(fid);
            filesWritten{end+1} = 'diagnostics_bysource.tsv';
        else
            warning('exportBayesianResults:OpenFailed', 'Could not open diagnostics_bysource.tsv');
        end
    end
end

%% 7. Posterior ensemble (accepted samples)
nAcc = 0;
if ~isempty(posteriorSamples) && isstruct(posteriorSamples)
    if isfield(posteriorSamples, 'kcats') && ~isempty(posteriorSamples.kcats)
        pk = posteriorSamples.kcats;
        nAcc = size(pk, 2);
        sampleHeaders = arrayfun(@(s) sprintf('sample%d',s), 1:nAcc, 'UniformOutput', false);
        pIDs = paramIDs; pSrc = paramSources;
        if size(pk,1) ~= D            % row labels only valid if dimension matches
            pIDs = {}; pSrc = {};
        end
        if writeParamMatrixTSV(fullfile(outputPath,'posterior_kcats.tsv'), pk, ...
                sampleHeaders, pIDs, pSrc, numFmtKcat)
            filesWritten{end+1} = 'posterior_kcats.tsv';
        end
    end
    if isfield(posteriorSamples, 'rmse') && ~isempty(posteriorSamples.rmse)
        pr = posteriorSamples.rmse(:);
        fid = fopen(fullfile(outputPath,'posterior_rmse.tsv'), 'w');
        if fid ~= -1
            fprintf(fid, 'sample_index\trmse\n');
            for s = 1:numel(pr)
                fprintf(fid, '%d\t%.8g\n', s, pr(s));
            end
            fclose(fid);
            filesWritten{end+1} = 'posterior_rmse.tsv';
        else
            warning('exportBayesianResults:OpenFailed', 'Could not open posterior_rmse.tsv');
        end
    end
end

%% 8. Run metadata (plain text key/value)
fid = fopen(fullfile(outputPath,'run_metadata.txt'), 'w');
if fid ~= -1
    fprintf(fid, 'bayesianSensitivityTuning export\n');
    fprintf(fid, 'generated\t%s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, 'n_parameters\t%d\n', D);
    fprintf(fid, 'n_trace_columns\t%d (generation 0..%d; column 1 = prior)\n', ...
        size(kcatTrace,2), size(kcatTrace,2)-1);
    fprintf(fid, 'n_generations_diagnostics\t%d\n', nDiag);
    if isfield(diagnostics, 'finalGeneration')
        fprintf(fid, 'final_generation\t%d\n', diagnostics.finalGeneration);
    end
    if isfield(diagnostics, 'converged')
        fprintf(fid, 'converged\t%d\n', double(diagnostics.converged));
    end
    if ~isempty(rmseVec)
        fprintf(fid, 'initial_rmse\t%.8g\n', rmseVec(1));
        fprintf(fid, 'final_rmse\t%.8g\n', rmseVec(end));
    end
    fprintf(fid, 'n_accepted_samples\t%d\n', nAcc);
    if isfield(diagnostics, 'sourceLabels')
        fprintf(fid, 'source_labels\t%s\n', strjoin(diagnostics.sourceLabels, ', '));
    end
    fprintf(fid, 'files_written\t%s\n', strjoin(filesWritten, ', '));
    fclose(fid);
    filesWritten{end+1} = 'run_metadata.txt';
end

fprintf('exportBayesianResults: wrote %d files to %s\n', numel(filesWritten), outputPath);

end


%% ------------------------------------------------------------------------
%% Local helpers
%% ------------------------------------------------------------------------

function ok = writeGenTable(filepath, genVec, headerCells, dataMat, numFmt)
% Generation-indexed table: column 1 = integer generation, then numeric cols.
% Returns true on success.
    ok = false;
    fid = fopen(filepath, 'w');
    if fid == -1
        warning('exportBayesianResults:OpenFailed', 'Could not open %s', filepath);
        return;
    end
    fprintf(fid, 'generation');
    for c = 1:numel(headerCells)
        fprintf(fid, '\t%s', headerCells{c});
    end
    fprintf(fid, '\n');
    for r = 1:numel(genVec)
        fprintf(fid, '%d', genVec(r));
        fprintf(fid, ['\t' numFmt], dataMat(r, :));
        fprintf(fid, '\n');
    end
    fclose(fid);
    ok = true;
end

function ok = writeParamMatrixTSV(filepath, dataMat, colHeaders, paramIDs, paramSources, numFmt)
% Parameter-indexed matrix: one row per parameter (model kcat entry).
%   Columns: param_index [, reaction_id] [, source], then colHeaders.
% Returns true on success.
    ok = false;
    [nRows, C] = size(dataMat);
    fid = fopen(filepath, 'w');
    if fid == -1
        warning('exportBayesianResults:OpenFailed', 'Could not open %s', filepath);
        return;
    end
    hasID  = ~isempty(paramIDs);
    hasSrc = ~isempty(paramSources);
    % Header
    fprintf(fid, 'param_index');
    if hasID,  fprintf(fid, '\treaction_id'); end
    if hasSrc, fprintf(fid, '\tsource');      end
    for c = 1:C
        fprintf(fid, '\t%s', colHeaders{c});
    end
    fprintf(fid, '\n');
    % Rows
    for r = 1:nRows
        fprintf(fid, '%d', r);
        if hasID,  fprintf(fid, '\t%s', paramIDs{r});     end
        if hasSrc, fprintf(fid, '\t%s', paramSources{r}); end
        fprintf(fid, ['\t' numFmt], dataMat(r, :));
        fprintf(fid, '\n');
    end
    fclose(fid);
    ok = true;
end