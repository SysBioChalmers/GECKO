function varargout = enzymeUsage(varargin)
% enzymeUsage  Deprecated alias for getEnzymeUsage.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to getEnzymeUsage.
%
% See also
% --------
% getEnzymeUsage

warning('GECKO:deprecatedName', ...
    'enzymeUsage is deprecated; use getEnzymeUsage instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = getEnzymeUsage(varargin{:});
end
