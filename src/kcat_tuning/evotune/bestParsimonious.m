function point = bestParsimonious(points,varargin)
% bestParsimonious  The fewest-changes point whose RMSE is within
% tolerance of the best on a parsimonyFrontier.
%
% Defaults to 2%, which is the scale of alternate-optimum noise typically
% seen from re-solving an LP with a slightly different kcat vector --
% below that, two vectors are not meaningfully different in fit, so the
% sparser one wins.
%
% Parameters
% ----------
% points : struct array
%     as returned by parsimonyFrontier.
%
% Name-Value Arguments
% --------------------
% tolerance : double
%     default 0.02.
%
% Returns
% -------
% point : struct
%     the chosen entry of points.
%
% See also
% --------
% parsimonyFrontier
p = parseGECKOargs(varargin, { ...
    'tolerance', 0.02});
tolerance = p.tolerance;

if isempty(points)
    error('points is empty.')
end
rmses = [points.rmse];
bestRmse = min(rmses);
within = points(rmses <= bestRmse * (1+tolerance));
[~,idx] = min([within.nChanged]);
point = within(idx);
end
