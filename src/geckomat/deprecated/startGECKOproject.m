function varargout = startGECKOproject(varargin)
% startGECKOproject  Deprecated alias for createGECKOproject.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to createGECKOproject.
%
% See also
% --------
% createGECKOproject

warning('GECKO:deprecatedName', ...
    'startGECKOproject is deprecated; use createGECKOproject instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = createGECKOproject(varargin{:});
end
