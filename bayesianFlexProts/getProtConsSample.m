function r = getProtConcSample(mu,sigma,step,method)
% getProtConcSample
%   Samples random protein concentrations from a distribution.
%
% Input:
%   mu              Mean of distribution (data is logged to get a normal distr)
%   sigma           Std deviation of the distribution
%   step            Number of protein concs to sample
%   method          shape of distribution: 'normal' or 'uniform'. 
%                   (Optional, default is 'normal')
% Output:
%   r               The sampled protein concs

if nargin < 4
    method = 'normal';
end
if mu == 0
    r = zeros(1,step);
elseif strcmp(method,'normal')
    mutmp = log10(mu);
    %sigmatmp = log10(sigma/3600);
    sigmatmp = sigma;
    pd = makedist('normal','mu',mutmp,'sigma',sigmatmp);
    %t = truncate(pd,-3,8);
    r = random(pd,1,step);
    r = 10.^(r);
elseif strcmp(method,'uniform')
    mutmp = log10(mu);
    sigmatmp = sigma;
    pd = makedist('uniform','lower',mutmp-sigmatmp,'upper',mutmp + sigmatmp);
    t = truncate(pd,-2,8);
    r = random(t,1,step);
    r = 10.^(r);
end

r(r<0) = 0;

end