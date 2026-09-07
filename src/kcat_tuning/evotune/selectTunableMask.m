function mask = selectTunableMask(ecModel,screen,varargin)
% selectTunableMask  Build a tunableMask from a screenKcatLeverage report.
%
% Keeps the fewest highest-ranked groups whose combined leverage reaches
% targetImpactShare of the total -- a relative cutoff, so the same
% targetImpactShare selects a comparable quality of parameter set on any
% model, rather than an absolute leverage value or a fixed count
% calibrated on a different model's scale.
%
% Pure and cheap -- no simulation, only screen is read.
%
% Parameters
% ----------
% ecModel : struct
% screen : table
%     as returned by screenKcatLeverage.
%
% Name-Value Arguments
% --------------------
% targetImpactShare : double
%     default 0.9. Not a validated universal constant -- check the
%     resulting mask's size against screenKcatLeverage's own curve before
%     trusting it on a new model.
%
% Returns
% -------
% mask : logical
%     length numel(ecModel.ec.rxns).
%
% See also
% --------
% screenKcatLeverage, cmaesKcatTuning
p = parseGECKOargs(varargin, { ...
    'targetImpactShare', 0.9});
targetImpactShare = p.targetImpactShare;

mask = false(numel(ecModel.ec.rxns),1);
if isempty(screen) || height(screen) == 0
    return
end
prevCum = [0; screen.cumLeverageShare(1:end-1)];
keep = prevCum < targetImpactShare;
for i = find(keep)'
    mask(screen.positions{i}) = true;
end
end
