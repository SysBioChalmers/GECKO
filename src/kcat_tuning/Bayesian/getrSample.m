function r = getrSample(mu,sigma,n,varargin)
% getrSample  Sample random kcat values from a distribution.
%
% Draws n random kcat values from a distribution defined by a mean and
% standard deviation. If mu is zero, a vector of zeros is returned.
%
% Parameters
% ----------
% mu : double
%     mean of the distribution.
% sigma : double
%     standard deviation of the distribution.
% n : double
%     number of kcat values to sample.
%
% Name-Value Arguments
% --------------------
% method : char
%     shape of the distribution, 'lognormal' or 'uniform' (default
%     'lognormal'). If 'lognormal' is selected, mu and sigma are
%     automatically converted into log-space.
%
% Returns
% -------
% r : double
%     sample of kcat values.
%
% Examples
% --------
%     % the optional argument may be given positionally or as a name-value pair:
%     r = getrSample(mu, sigma, n);
%     r = getrSample(mu, sigma, n, 'method', 'uniform');
%     r = getrSample(mu, sigma, n, 'uniform');

p = parseGECKOargs(varargin, { ...
    'method', 'lognormal'});
method = p.method;

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
