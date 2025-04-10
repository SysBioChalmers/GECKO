function [a,b] = updatepriorProts(x)
% updatepriorProts
%   Calculates a new distribution from the selected protein concentration values
%
% Input:
%   x               concs
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
    pd = fitdist(log10(x),'Normal');
    a = 10^(pd.mu);
    b = pd.sigma;
end

end