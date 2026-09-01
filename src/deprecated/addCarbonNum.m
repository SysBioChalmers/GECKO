function varargout = addCarbonNum(varargin)
% addCarbonNum  Deprecated alias for fillCarbonNum.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to fillCarbonNum.
%
% See also
% --------
% fillCarbonNum

warning('GECKO:deprecatedName', ...
    'addCarbonNum is deprecated; use fillCarbonNum instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = fillCarbonNum(varargin{:});
end
