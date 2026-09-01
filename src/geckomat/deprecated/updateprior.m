function varargout = updateprior(varargin)
% updateprior  Deprecated alias for updatePrior.
%
% Kept for backward compatibility with pre-GECKO4 code. Will be removed
% in a future release; switch to updatePrior.
%
% See also
% --------
% updatePrior

warning('GECKO:deprecatedName', ...
    'updateprior is deprecated; use updatePrior instead. The old name will be removed in a future release.');
[varargout{1:nargout}] = updatePrior(varargin{:});
end
