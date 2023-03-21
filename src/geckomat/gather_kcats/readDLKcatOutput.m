function kcatList = readDLKcatOutput(model, outFile, modelAdapter)
% readDLKcatOutput
%   Reads the DLKcat output file and constructs a kcatList structure, that
%   can be used by selectKcatValue() to populate the ecModel with kcat
%   values.
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   outFile         name and path of the DLKcat output file. If nothing is
%                   provided, an attempt will be made to read
%                   data/DLKcat.tsv from the obj.params.path folder
%                   specified in the modelAdapter.
%   modelAdapter    a loaded model adapter (Optional, will otherwise use the
%                   default model adapter).
%
% Output:
%   kcatList    structure array with list of DLKcat derived kcat values,
%               with separate entries for each kcat value
%               source      'DLKcat'           
%               rxns        reaction identifiers
%               genes       gene identifiers
%               substrate   substrate names
%               kcat        predicted kcat value in /sec
%
% Usage:
%   kcatList = readDLKcatOutput(model, outFile, modelAdapter)

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if nargin<2 || isempty(outFile)
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

% Check that all substrates are in the model
if ~all(ismember(subs,model.metNames))
    error('Not all substrates from DLKcat output can be found in model.metNames. DLKcat was likely run with an input file that was generated from another ecModel.')
end

% Check that all reactions are in model.ec.rxns
if ~all(ismember(rxns,model.ec.rxns))
    error('Not all reactions from DLKcat output can be found in model.ec.rxns. DLKcat was likely run with an input file that was generated from another ecModel.')
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
