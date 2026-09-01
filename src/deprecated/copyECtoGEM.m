function varargout = copyECtoGEM(varargin)
% copyECtoGEM  Deprecated alias for applyECcodes.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to applyECcodes.
%
% See also
% --------
% applyECcodes

warning('GECKO:deprecatedName', ...
    'copyECtoGEM is deprecated; use applyECcodes instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = applyECcodes(varargin{:});
end
