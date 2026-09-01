function varargout = selectKcatValue(varargin)
% selectKcatValue  Deprecated alias for assignKcatValues.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to assignKcatValues.
%
% See also
% --------
% assignKcatValues

warning('GECKO:deprecatedName', ...
    'selectKcatValue is deprecated; use assignKcatValues instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = assignKcatValues(varargin{:});
end
