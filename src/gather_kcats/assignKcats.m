function [model, decisions, report] = assignKcats(model, varargin)
% assignKcats  Assign kcats from the kcat source files and apply them.
%
% Decides one kcat per reaction from the kcat source files in data/kcats,
% as [kcat_assignment] in model_adapter.toml sets (see collectKcats), writes
% model.ec.kcat and model.ec.source, and applies the kcat constraints
% (applyKcatConstraints). Any kcats the model had are replaced: the source
% files are the input, and nothing in the model is written back to them.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% writeReport : logical
%     write the decisions to output/<model id>_kcatDecisions.tsv (default
%     true).
%
% Returns
% -------
% model : struct
%     the ecModel with the assigned kcats applied.
% decisions, report
%     as returned by collectKcats.
%
% Examples
% --------
%     updateKcatSources(ecModel);   % gather what the configuration needs
%     ecModel = assignKcats(ecModel);
%
% See also
% --------
% collectKcats, updateKcatSources, applyKcatConstraints

checkNameValue(varargin, {'modelAdapter','writeReport'}, 'assignKcats');
[decisions, report, model] = collectKcats(model, varargin{:});
model = applyKcatConstraints(model);
end
