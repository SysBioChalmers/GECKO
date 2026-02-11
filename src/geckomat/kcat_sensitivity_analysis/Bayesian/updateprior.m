function [mu,sigma] = updateprior(x,defaultCV)
% updateprior
%   Calculates a new distribution from the provided kcat values
%
% Input:
%   x       kcat values
%
% Output:
%   mu      mean
%   sigma   standard deviation

if nargin<2
    defaultCV = 0.25;
end

if iscell(x), x = x{:}; end

% Remove non-positive values—they cannot be logged
n_before = numel(x);
x = x(x > 0);
if isempty(x)
    error('x cannot be empty or contain only non-positive values.');
end
if numel(x) < n_before
    warning('Removed %d non-positive values before fitting.', n_before - numel(x));
end

if isscalar(x)
    mu      = x;
    sigma   = x * defaultCV;
else
    pd      = fitdist(log(x),'normal');
    mu    = exp(pd.mu + 0.5 * pd.sigma^2);
    sigma = sqrt( (exp(pd.sigma^2) - 1) * exp(2*pd.mu + pd.sigma^2) );
end
end
