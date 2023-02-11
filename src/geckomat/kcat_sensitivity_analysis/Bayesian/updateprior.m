function [a,b] = updateprior(x)
% updateprior
%   Calculates a new distribution from the selected kcat values
%
% Input:
%   x               kcats
% Output:
%   a               Mean
%   b               Std dev

x = x{:};
x(x == 0) = []; %remove zeros - they cannot be handled in the log transform below
if length(x) == 0
    %we don't have much choice in the two first cases - just set sigma to 1 - same as in the initial prior
    a = 0;
    b = 1;
elseif length(x) == 1
    a = x;
    b = 1;
else
    pd = fitdist(log10(x./3600),'Normal');
    a = 10^(pd.mu)*3600;
    b = pd.sigma;
end

end