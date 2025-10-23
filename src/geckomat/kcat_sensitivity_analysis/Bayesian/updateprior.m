function [mu,sigma] = updateprior(x)
% updateprior
%   Calculates a new distribution from the provided kcat values
%
% Input:
%   x       kcat values
%
% Output:
%   mu      mean
%   sigma   standard deviation

x = x{:};
%x(x == 0) = []; %remove zeros - they cannot be handled in the log transform below
if isempty(x)
    error('x cannot be empty')
elseif isscalar(x)
    % By default set sigma at 25% of mu
    mu      = x;
    sigma   = x/4;
else
    pd      = fitdist(log(x),'normal');
    mu      = exp(pd.mu);
    sigma   = mu * sqrt(exp(pd.sigma^2) - 1);
end
end
