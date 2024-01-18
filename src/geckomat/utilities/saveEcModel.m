function saveEcModel(ecModel,filename,modelAdapter)
% saveECmodel
%   Saves the ecModel in either YAML format (= default and preferred, as
%   all ecModel content is reserved) and/or SBML format (= more widely
%   compatible with other constraint-based modelling tools, can be used for
%   running simulations like FBA etc., but this model cannot be loaded back
%   into MATLAB for applying further GECKO functions, as model.ec is lost).
%
%   The ecModel is saved to the ecModel-specific models/ subfolder. For
%   saving to other locations, writeYAMLmodel or exportModel are more
%   suitable.
%
% Input:
%   ecModel         an ecModel in GECKO 3 format (with ecModel.ec structure)
%   filename        ending with either .yml or .xml, specifying if the
%                   ecModel should be saved in YAML or SBML file format. If
%                   no file extension is given, the ecModel will be saved
%                   in YAML format. If no filename is given, 'ecModel.yml'
%                   is used.
%   modelAdapter    a loaded model adapter, from where the model folder is
%                   read (Optional, will otherwise use the default model adapter).
%
% Usage:
%   saveECmodel(ecModel,filename,modelAdapter)


if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
if nargin < 2 || isempty(filename)
    filename = 'ecModel';
end
filename = fullfile(params.path,'models',filename);

ecModel.description = ['Enzyme-constrained model of ' ecModel.id];

switch filename(end-3:end)
    case {'.xml','sbml'}
        exportModel(ecModel, filename);
        %Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
        %inconsistencies between Windows and MAC:
        copyfile(filename,'backup.xml')
        fin  = fopen('backup.xml', 'r');
        fout = fopen(filename, 'w');
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
    otherwise % Assume YAML
        writeYAMLmodel(ecModel,filename);
end
end
