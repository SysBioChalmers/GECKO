function checkNameValue(args, names, functionName)
% checkNameValue  Accept only name-value pairs with known names.
%
% The kcat source functions have no positional form, so a misspelled
% option name must be an error rather than being read as a positional
% argument by parseGECKOargs.
id = [functionName ':arguments'];
if mod(numel(args), 2) ~= 0
    error(id, 'Options of %s are given as name-value pairs.', functionName)
end
for k = 1:2:numel(args)
    name = args{k};
    if ~(ischar(name) || (isstring(name) && isscalar(name)))
        error(id, 'Option names of %s are text, not %s. Options: %s.', functionName, ...
            class(name), strjoin(names, ', '))
    elseif ~any(strcmp(char(name), names))
        error(id, 'Unknown option of %s: %s. Options: %s.', functionName, char(name), ...
            strjoin(names, ', '))
    end
end
end
