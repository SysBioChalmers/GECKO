function r = getrSample(mu,sigma,step,method)
% getrSample
%   Samples random kcats from a distribution.
%
% Input:
%   mu              Mean of distribution (data is logged to get a normal distr)
%   sigma           Std deviation of the distribution
%   step            Number of kcats to sample
%   method          shape of distribution: 'normal' or 'uniform'. 
%                   (Optional, default is 'normal')
% Output:
%   r               The sampled kcats
%
if nargin < 4
    method = 'normal';
end
if mu == 0
    r = zeros(1,step);
elseif strcmp(method,'normal')
    mutmp = log10(mu/3600);
    %sigmatmp = log10(sigma/3600);
    sigmatmp = sigma;
    pd = makedist('normal','mu',mutmp,'sigma',sigmatmp);
    %t = truncate(pd,-3,8);
    r = random(pd,1,step);
    r = 10.^(r).*3600;
elseif strcmp(method,'uniform')
    mutmp = log10(mu/3600);
    sigmatmp = sigma;
    pd = makedist('uniform','lower',mutmp-sigmatmp,'upper',mutmp + sigmatmp);
    t = truncate(pd,-2,8);
    r = random(t,1,step);
    r = 10.^(r).*3600;
end

r(r<0) = 0;

end