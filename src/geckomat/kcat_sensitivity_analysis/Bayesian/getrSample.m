function r = getrSample(mu,sigma,n,method)
% getrSample
%   Sample random kcat values from a distribution.
%
% Input:
%   mu      mean of distribution
%   sigma   standard deviation of distribution
%   n       number of kcat values to sample
%   method  shape of distribution: 'lognormal' or 'uniform'. (optional,
%           default is 'lognormal'). If 'lognormal' is selected, mu and
%           sigma are automatically converted into log-space
%
% Output:
%   r       sample of kcat values

if nargin < 4
    method = 'lognormal';
end
if mu == 0
    r = zeros(1,n);
elseif strcmp(method,'lognormal')
    % Convert to log-space
    varLog      = log(1 + (sigma/mu)^2);
    sigmaLog    = sqrt(varLog);
    muLog       = log(mu) - 0.5 * varLog;

    % Take random samples from lognormal distribution
    r = lognrnd(muLog,sigmaLog,[1,n]);
elseif strcmp(method,'uniform')
    pd = makedist('uniform','lower',mu - sigma,'upper',mu + sigma);
    t = truncate(pd,0,inf);
    r = random(t,1,n);
end
end
