function args = parseGECKOargs(rawArgs, spec)
% parseGECKOargs  Resolve optional arguments passed positionally or by name.
%
% Helper that lets a GECKO function accept its optional arguments either in
% the traditional positional order or as name-value pairs, without losing
% backward compatibility. The required leading arguments stay explicit in the
% function signature; the optional tail is collected into varargin and passed
% here together with a specification of the optional parameters.
%
% Parameters
% ----------
% rawArgs : cell
%     the function's varargin, i.e. whatever was supplied after the required
%     positional arguments.
% spec : cell
%     an N-by-2 or N-by-3 cell array describing the optional parameters, one
%     per row:
%
%     - column 1 : parameter name (char).
%     - column 2 : default value, used when the parameter is not supplied.
%     - column 3 : optional validator, a function handle that is called with
%       the resolved value and should error on an invalid value. Use [] for
%       no validation.
%
% Returns
% -------
% args : struct
%     a struct with one field per parameter name, holding the resolved value.
%
% Notes
% -----
% The calling convention is detected from rawArgs: it is treated as name-
% value pairs only when it has an even number of elements AND every element
% in a name position is a string matching one of the declared parameter
% names. Otherwise it is treated as positional and assigned in the order of
% the spec. A call therefore uses either the positional or the name-value
% form for the optional tail, not a mix of the two. When an optional tail
% value could itself be a string equal to a parameter name, use the explicit
% name-value form to remove the ambiguity.
%
% Examples
% --------
%     % inside makeEcModel(model, varargin):
%     p = parseGECKOargs(varargin, { ...
%         'geckoLight',   false; ...
%         'modelAdapter', []});
%     geckoLight   = p.geckoLight;
%     modelAdapter = p.modelAdapter;
%
%     % the caller may then use either form:
%     makeEcModel(model, true, adapter)                  % positional
%     makeEcModel(model, 'geckoLight', true)             % named

if nargin < 2 || isempty(spec)
    error('parseGECKOargs:noSpec', 'A parameter specification is required.');
end

names = spec(:, 1);
defaults = spec(:, 2);
hasValidators = size(spec, 2) >= 3;

% Seed every parameter with its default.
args = struct();
for i = 1:numel(names)
    args.(names{i}) = defaults{i};
end

if ~isempty(rawArgs)
    isNameValue = mod(numel(rawArgs), 2) == 0 && ...
        all(cellfun(@(x) (ischar(x) || isstring(x)) && isscalar(string(x)) && ...
        any(strcmpi(x, names)), rawArgs(1:2:end)));

    if isNameValue
        for k = 1:2:numel(rawArgs)
            args.(names{strcmpi(rawArgs{k}, names)}) = rawArgs{k + 1};
        end
    else
        if numel(rawArgs) > numel(names)
            error('parseGECKOargs:tooManyArgs', ...
                'Too many positional arguments: expected at most %d, got %d.', ...
                numel(names), numel(rawArgs));
        end
        for k = 1:numel(rawArgs)
            args.(names{k}) = rawArgs{k};
        end
    end
end

% Apply optional validators.
if hasValidators
    for i = 1:numel(names)
        validator = spec{i, 3};
        if ~isempty(validator)
            validator(args.(names{i}));
        end
    end
end
end
