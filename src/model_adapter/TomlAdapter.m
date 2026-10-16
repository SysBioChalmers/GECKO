% TomlAdapter  Model adapter configured by model_adapter.toml.
%
% The model adapter for a project whose parameters are in
% model_adapter.toml, the configuration file that GECKO and geckopy share.
% No adapter class has to be written: what used to be code in an adapter
% class is set in the TOML file.
%
% - Spontaneous reactions: [spontaneous] name_contains (reactions whose
%   name contains the text) and/or file (a TSV of reaction identifiers,
%   optionally with a second 0/1 column).
% - Anaerobic conditions and protein content during kcat tuning:
%   [matlab] make_model_anaerobic and change_protein_biomass name MATLAB
%   functions, for instance in the project's code folder, which is added
%   to the MATLAB path. Each takes the ecModel (and, for
%   change_protein_biomass, the protein content) and returns the ecModel.
%
% A project that needs other custom code can keep an adapter class and
% read its parameters from the TOML file in its constructor:
%
%     obj.params = loadAdapterParams(fullfile(projectPath,'model_adapter.toml'));
%
% Examples
% --------
%     adapter = TomlAdapter(fullfile(projectPath,'model_adapter.toml'));
%     % or, as the default adapter:
%     ModelAdapterManager.setDefault(projectPath);
%
% See also
% --------
% ModelAdapter, ModelAdapterManager, loadAdapterParams, convertAdapterToToml
classdef TomlAdapter < ModelAdapter
    methods
        function obj = TomlAdapter(tomlFile)
            % TomlAdapter  Load the adapter from model_adapter.toml.
            %
            % Parameters
            % ----------
            % tomlFile : char
            %     path to model_adapter.toml, or to the project folder that
            %     holds it.
            obj.params = loadAdapterParams(tomlFile);
            codeFolder = fullfile(obj.params.path, 'code');
            if isfolder(codeFolder)
                addpath(codeFolder);
            end
        end

        function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
            % getSpontaneousReactions  Reactions set as spontaneous in [spontaneous].
            %
            % Returns
            % -------
            % spont : logical
            %     true for each reaction of model that is spontaneous.
            % spontRxnNames : cell
            %     identifiers of those reactions.
            spont = false(numel(model.rxns),1);
            if isfield(obj.params,'spontaneous')
                s = obj.params.spontaneous;
                if isfield(s,'name_contains') && ~isempty(s.name_contains)
                    spont = spont | contains(model.rxnNames, s.name_contains);
                end
                if isfield(s,'file') && ~isempty(s.file)
                    spont = spont | ismember(model.rxns, readSpontaneousFile(s.file));
                end
            end
            spontRxnNames = model.rxns(spont);
        end

        function ecModel = makeModelAnaerobic(obj,ecModel)
            % makeModelAnaerobic  Run the function named in [matlab] make_model_anaerobic.
            fun = hookName(obj, 'makeModelAnaerobic');
            if ~isempty(fun)
                ecModel = feval(fun, ecModel);
            end
        end

        function ecModel = changeProteinBiomass(obj,ecModel,Ptot)
            % changeProteinBiomass  Run the function named in [matlab] change_protein_biomass.
            fun = hookName(obj, 'changeProteinBiomass');
            if ~isempty(fun)
                ecModel = feval(fun, ecModel, Ptot);
            end
        end
    end
end

function fun = hookName(obj, field)
fun = '';
if isfield(obj.params,'matlab') && isfield(obj.params.matlab, field)
    fun = obj.params.matlab.(field);
end
end

function rxns = readSpontaneousFile(file)
% Reaction identifiers in the first column; with a second column, only
% rows where it is 1 or true. The first row is a header.
fid = fopen(file, 'r');
if fid < 0
    error('GECKO:TomlAdapter', 'Cannot read the spontaneous reactions file %s', file);
end
header = fgetl(fid);
nCols = numel(strsplit(header, '\t'));
if nCols >= 2
    cols = textscan(fid, ['%q %q' repmat(' %*q', 1, nCols-2)], 'Delimiter', '\t');
    keep = ismember(lower(strtrim(cols{2})), {'1','true'});
    rxns = cols{1}(keep);
else
    cols = textscan(fid, '%q', 'Delimiter', '\t');
    rxns = cols{1};
end
fclose(fid);
end
