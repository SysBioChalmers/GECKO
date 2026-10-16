function kcatList = kcatsFromDLKcat(model, varargin)
% kcatsFromDLKcat  Gather kcats from DLKcat and keep them as a source file.
%
% Predicts kcats with DLKcat and writes them to data/kcats/dlkcat.tsv in
% the project folder (see saveKcatSource). With 'run', true it writes the
% DLKcat input (writeDLKcatInput) and runs DLKcat through its Docker image
% (runDLKcat), which needs Docker; otherwise it reads the output of an
% earlier run (readDLKcatOutput). A run records how many reactions it
% covered in data/DLKcat_run.txt, so that the source file says so also when
% the output is read later.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure), with
%     metabolite SMILES (findMetSmiles) and protein sequences.
%
% Name-Value Arguments
% --------------------
% run : logical
%     write the input and run DLKcat; false reads the existing output in
%     data/DLKcat.tsv (default false).
% ecRxns : logical
%     for which model.ec.rxns DLKcat is run (default: all). Only used with
%     'run', true.
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'dlkcat').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values, as returned by readDLKcatOutput.
%
% Examples
% --------
%     kcatsFromDLKcat(ecModel, 'run', true);
%     kcatList = kcatsFromDLKcat(ecModel);
%
% See also
% --------
% writeDLKcatInput, runDLKcat, readDLKcatOutput, kcatsFromBrenda, saveKcatSource

checkNameValue(varargin, {'run','ecRxns','save','saveAs','modelAdapter'}, 'kcatsFromDLKcat');
p = parseGECKOargs(varargin, { ...
    'run',          false; ...
    'ecRxns',       []; ...
    'save',         true; ...
    'saveAs',       'dlkcat'; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
if ~p.run && ~isempty(p.ecRxns)
    error('kcatsFromDLKcat:ecRxns', ['ecRxns selects the reactions DLKcat is run for, so ' ...
        'it needs ''run'', true.'])
end
if p.save
    file = kcatSourceFile(p.saveAs, [], modelAdapter);
    if p.run
        warnSubset(p.ecRxns, numel(model.ec.rxns), file);
    end
end

runFile = fullfile(modelAdapter.params.path, 'data', 'DLKcat_run.txt');
if p.run
    % As runDLKcat: on macOS, MATLAB started from the launcher may not
    % have Docker on its path.
    if ismac
        setenv('PATH', strcat('/usr/local/bin', ':', getenv("PATH")));
    end
    [status, ~] = system('docker --version');
    if status ~= 0
        error('kcatsFromDLKcat:noDocker', ['DLKcat runs in a Docker image, but Docker was ' ...
            'not found. Install it from https://www.docker.com, or call kcatsFromDLKcat ' ...
            'without ''run'' to read an earlier DLKcat output.'])
    end
    writeDLKcatInput(model, 'ecRxns', p.ecRxns, 'modelAdapter', modelAdapter, 'overwrite', true);
    runDLKcat('modelAdapter', modelAdapter);
    fid = fopen(runFile, 'w');
    if ~isempty(p.ecRxns)
        fprintf(fid, 'ecRxns: %d\n', countRxns(p.ecRxns));
    end
    fclose(fid);
else
    outFile = fullfile(modelAdapter.params.path, 'data', 'DLKcat.tsv');
    if ~isfile(outFile)
        error('kcatsFromDLKcat:noOutput', ['No DLKcat output at %s. Run ' ...
            'kcatsFromDLKcat(model, ''run'', true) first.'], outFile)
    end
end
kcatList = readDLKcatOutput(model, 'modelAdapter', modelAdapter);

if p.save
    parameters = struct('run', logical(p.run));
    if isfile(runFile)
        tok = regexp(fileread(runFile), 'ecRxns:\s*(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            parameters.ec_rxns = str2double(tok{1});
        end
    end
    saveKcatSource(model, kcatList, p.saveAs, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromDLKcat', 'parameters', parameters);
    fprintf('DLKcat kcats for %d reactions written to %s\n', ...
        numel(unique(kcatList.rxns(kcatList.kcats > 0))), file);
end
end
