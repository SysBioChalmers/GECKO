function model = loadEcModel(filename, modelAdapter)
% loadEcModel
%   Loads the ecModel that matches the modelAdapter. By default, it loads
%   the models/ecModel.yml in the directory specified as param.path in
%   the modelAdapter. Alternative files in the same folder can be loaded by
%   providing the appropriate filename. If loading models from other
%   locations, one can directly use readYAMLmodel.
%
% Input:
%   filename        name of the ecModel file (Optional, default 'ecModel.yml').
%                   Fill should be located in the models/ subfolder of
%                   param.path as specified in the modelAdapter.
%   modelAdapter    a loaded model adapter, from where the model folder is
%                   read (Optional, will otherwise use the default model adapter).
%
% Output:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%
% Usage:
%   model = loadEcModel(filename, modelAdapter)

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
if nargin < 1 || isempty(filename)
    filename = 'ecModel.yml';
else
end

if ~endsWith(filename,{'yml','yaml'})
    error(['ecModels should preferably be distributed in YAML file format, ', ...
           'as otherwise model content will be lost. If the model is in ', ...
           'SBML format, you can use the usual importModel function to load ', ...
           'this as any other genome-scale model.'])
end

filename = fullfile(params.path,'models',filename);
if endsWith(filename,{'yml','yaml'})
    model = readYAMLmodel(filename);
elseif endsWith(filename,{'xml','sbml'})
    model = importModel(filename);
end

% Older ecModels (built before this version of GECKO) had their
% protein usage and protein pool reactions defined as "reverse"
% reactions: their flux was negative, with the usual stoichiometry
% signs swapped accordingly. Newer GECKO builds them as ordinary
% forward reactions (positive flux). Detect the older convention
% when loading and flip those reactions in place, so the rest of
% the GECKO functions always see the same shape.
model = flipLegacyProtDirection(model);
end

function model = flipLegacyProtDirection(model)
%Detect whether the protein-related reactions in `model` follow the
%older "reverse direction" convention, and flip them to the current
%forward convention if so. The check looks for any usage_prot_* or
%prot_pool_exchange reaction whose lower bound is negative.
usageIdx = find(startsWith(model.rxns,'usage_prot_'));
poolIdx  = find(strcmp(model.rxns,'prot_pool_exchange'));
flipIdx  = [usageIdx; poolIdx];
flipIdx  = flipIdx(model.lb(flipIdx) < -1e-9);
if isempty(flipIdx)
    return
end
printOrange(['INFO: ecModel uses the older reverse-direction convention ' ...
             'for protein usage / pool reactions. Flipping them to the ' ...
             'current forward convention.\n']);
%Flip the stoichiometry signs and swap the bounds for the affected
%reactions, so what was a reverse reaction becomes the equivalent
%forward one.
model.S(:,flipIdx) = -model.S(:,flipIdx);
oldLb = model.lb(flipIdx);
oldUb = model.ub(flipIdx);
model.lb(flipIdx) = -oldUb;
model.ub(flipIdx) = -oldLb;
if isfield(model,'rev')
    %Flux-only-forward reactions (lb >= 0) get rev=0.
    model.rev(flipIdx(model.lb(flipIdx) >= -1e-9)) = 0;
end
end
