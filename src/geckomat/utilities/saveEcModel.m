function saveECmodel(model,modelAdapter,format,filename)
% saveECmodel
%   Saves the model in either YAML format (= default and preferred, as all
%   model content is reserved) and/or SBML format (= more widely compatible
%   with other constraint-based modelling tools, can be used for running
%   simulations like FBA etc., but this model cannot be loaded back into
%   MATLAB for applying further GECKO functions, as model.ec is lost).
%
% Input:
%   model           an ecModel in GECKO 3 format (with ecModel.ec structure)
%   modelAdapter    a loaded model adapter, from where the model folder is
%                   read (Optional, will otherwise use the default model adapter).
%   format          file format, either as string if the model should be
%                   saved in one format, e.g. 'yaml', or cell array if
%                   multiple formats {'yaml','sbml'}. (Optional, default
%                   'yaml').
%   filename        overwrites whatever path is specified in modelAdapter,
%                   or if no modelAdapter is specified. Otherwise, the file
%                   is stored as ecModel in param.path.
%
% Usage:
%   saveECmodel(model,modelAdapter,format,filename)


if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
if nargin < 3 || isempty(format)
    format = 'yaml';
end
if nargin < 4 || isempty (format)
    filename = 'ecModel';
end
filename = fullfile(params.path,'models',filename);

%Model description:
model.description = ['Enzyme-constrained model of ' model.id];

if contains(format,{'yaml','yml'})
    writeYAMLmodel(model,[filename '.yml']);
end

if contains(format,{'xml','sbml'})
    exportModel(model,[filename '.xml']);
    %Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
    %inconsistencies between Windows and MAC:
    copyfile([filename '.xml'],'backup.xml')
    fin  = fopen('backup.xml', 'r');
    fout = fopen([filename '.xml'], 'w');
    still_reading = true;
    while still_reading
        inline = fgets(fin);
        if ~ischar(inline)
            still_reading = false;
        else
            if ~isempty(regexp(inline,'-00[0-9]','once'))
                inline = strrep(inline,'-00','-0');
            elseif ~isempty(regexp(inline,'-01[0-9]','once'))
                inline = strrep(inline,'-01','-1');
            end
            fwrite(fout, inline);
        end
    end
    fclose('all');
    delete('backup.xml');
end
end
