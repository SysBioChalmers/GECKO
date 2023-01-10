function kcatList = readDLKcatOutput(model,outFile)
% readDLKcatOutput
%   Reads the DLKcat output file and constructs a kcatList structure, that
%   can be used by selectKcatValue() to populate the ec-model with kcat
%   values.
%
% Input:
%   model       an ec-model in RAVEN format
%   outFile     name and path of the DLKcat output file. If nothing is
%               provided, an attempt will be made to read
%               data/DLKcatOutput.tsv from the obj.params.path folder
%               specified in the ModelAdapter.
%
% Output:
%   kcatList    structure array with list of DLKcat derived kcat values,
%               with separate entries for each kcat value
%               source      'DLKcat'           
%               rxns        reaction identifiers
%               genes       gene identifiers
%               substrate   substrate names
%               kcat        predicted kcat value in /sec

global GECKOModelAdapter
params=checkGECKOModelAdapter(GECKOModelAdapter);

if nargin<2
    fID      = fopen(fullfile(params.path,'data','DLKcatOutput.tsv'),'r');
else
    fID      = fopen(outFile);
end
DLKcatOutput = textscan(fID,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fID);

% Check that DLKcat output file and model match (not fool proof, but good enough)
[rxns, genes, subs, kcats] = deal(DLKcatOutput{[1,2,3,6]});

% Check that all substrates are in the model
if ~all(ismember(subs,model.metNames))
    error('Not all substrates from DLKcat output can be found in model.metNames')
end
% Check that all reactions are in model.ec.rxns
if ~all(ismember(rxns,model.ec.rxns))
    error('Not all reactions from DLKcat output can be found in model.ec.rxns')
end

% Filter out entries with no kcat value
noOutput        = strcmp(kcats,'None');
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
