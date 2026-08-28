function kcatList = readDLKcatOutput(model, varargin)
% readDLKcatOutput  Read the DLKcat output file into a kcatList structure.
%
% Reads the DLKcat output file and constructs a kcatList structure, that can
% be used by selectKcatValue() to populate the ecModel with kcat values.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% outFile : char
%     name and path of the DLKcat output file (default: data/DLKcat.tsv from
%     the obj.params.path folder specified in the modelAdapter).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     structure array with list of DLKcat derived kcat values, with separate
%     entries for each kcat value.
%
% Notes
% -----
% The kcatList structure has the following fields:
%
% - source : 'DLKcat'.
% - rxns : reaction identifiers.
% - genes : gene identifiers.
% - substrate : substrate names.
% - kcat : predicted kcat value in /sec.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     kcatList = readDLKcatOutput(model);
%     kcatList = readDLKcatOutput(model, 'outFile', 'myDLKcat.tsv');
%
% See also
% --------
% runDLKcat, writeDLKcatInput, selectKcatValue

p = parseGECKOargs(varargin, { ...
    'outFile',      []; ...
    'modelAdapter', []});
outFile      = p.outFile;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(outFile)
    fID      = fopen(fullfile(params.path,'data','DLKcat.tsv'),'r');
else
    fID      = fopen(outFile);
end
DLKcatOutput = textscan(fID,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fID);

% Check that DLKcat output file and model match (not fool proof, but good enough)
[rxns, genes, subs, kcats] = deal(DLKcatOutput{[1,2,3,6]});

% Check if it contains any kcat values
if all(cellfun(@isempty,kcats)) || all(strcmpi(kcats,'NA'))
    error('DLKcat file does not contain any kcat values, please run runDLKcat() first.')
end

% Check that all substrates are in the model (case-insensitive: an SBML loader can
% normalise metabolite-name capitalisation differently than whatever produced the
% DLKcat input, and that must not read as a genuinely missing substrate).
matchMets = ismember(lower(subs),lower(model.metNames));
if ~all(matchMets)
    errorText = 'DLKcat was likely run with an input file that was generated from another ecModel, as the following substrates from DLKcat output cannot be found in model.metNames:';
    error('%s', ravenList(errorText,subs(~matchMets),false))
end

% Check that all reactions are in model.ec.rxns
matchRxns = ismember(rxns,model.ec.rxns);
if ~all(matchRxns)
    errorText = 'DLKcat was likely run with an input file that was generated from another ecModel, as the following reactions from DLKcat output cannot be found in model.metNames:';
    error('%s', ravenList(errorText,rxns(~matchRxns),false))
end

% Filter out entries with no numeric value
noOutput        = cellfun(@isempty,regexpi(kcats,'[0-9]'));
kcats           = str2double(kcats(~noOutput));
rxns(noOutput)  = [];
genes(noOutput) = [];
subs(noOutput)  = [];

% Make kcatList structure
kcatList.source     = 'DLKcat';
kcatList.rxns       = rxns;
kcatList.genes      = genes;
kcatList.substrates = subs;
kcatList.kcats      = kcats;
end
