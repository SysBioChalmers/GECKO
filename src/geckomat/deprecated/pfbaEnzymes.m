function varargout = pfbaEnzymes(varargin)
% pfbaEnzymes  Deprecated alias for getPfbaEnzymes.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to getPfbaEnzymes.
%
% See also
% --------
% getPfbaEnzymes

warning('GECKO:deprecatedName', ...
    'pfbaEnzymes is deprecated; use getPfbaEnzymes instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = getPfbaEnzymes(varargin{:});
end
