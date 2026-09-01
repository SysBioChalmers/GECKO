function varargout = updateGECKOdoc(varargin)
% updateGECKOdoc  Deprecated alias for buildGECKOdoc.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to buildGECKOdoc.
%
% See also
% --------
% buildGECKOdoc

warning('GECKO:deprecatedName', ...
    'updateGECKOdoc is deprecated; use buildGECKOdoc instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = buildGECKOdoc(varargin{:});
end
