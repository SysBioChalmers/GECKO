function n = nChanged(kcat,kcat0,varargin)
% nChanged  How many kcats differ from their prior by more than relTol.
%
% The 2% default matches the threshold MATLAB's own per-source
% "unchanged" reporting has historically used.
%
% Parameters
% ----------
% kcat, kcat0 : double
%     tuned and prior kcat vectors, same length.
%
% Name-Value Arguments
% --------------------
% relTol : double
%     default 0.02.
%
% Returns
% -------
% n : double
%     count of entries where abs(log(kcat./kcat0)) > log1p(relTol).
%
% See also
% --------
% tunePriorPenaltyWeight, parsimonyFrontier
p = parseGECKOargs(varargin, { ...
    'relTol', 0.02});
relTol = p.relTol;
n = sum(abs(log(kcat(:) ./ kcat0(:))) > log1p(relTol));
end
