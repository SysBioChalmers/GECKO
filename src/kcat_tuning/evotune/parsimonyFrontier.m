function points = parsimonyFrontier(scoreFcn,kcat,kcat0,sigma0log,varargin)
% parsimonyFrontier  Score kcat with progressively more of it reverted to
% prior.
%
% A tuned kcat vector that fits well by rewriting every parameter is not
% a useful answer: changes should be few, and concentrated in the
% parameters least confident about. This sweeps a revert-threshold (in
% units of each kcat's own prior sigma) and reports how RMSE trades off
% against how many kcats stayed at their tuned value, so a sparser
% result within noise of the best fit can be preferred over the raw
% tuner output. See bestParsimonious to pick a point off this frontier.
%
% Parameters
% ----------
% scoreFcn : function_handle
%     scores one kcat vector; called once per threshold. Build it against
%     a specific ecModel/evoData, e.g.
%         scoreFcn = @(vec) evotuneScoreVector(ecModel,tunableIdx, ...
%             ecRxnIdsTunable,evoData,modelAdapter,params,vec);
%     where evotuneScoreVector applies vec to ecModel.ec.kcat(tunableIdx),
%     calls applyKcatConstraints, and returns evotuneScore's plain rmse.
% kcat, kcat0, sigma0log : double
%     tuned vector, prior vector, and per-kcat prior sigma, all the same
%     length.
%
% Name-Value Arguments
% --------------------
% thresholds : double
%     sigma cutoffs to sweep, ascending (default
%     [0 1 2 3 4 6 8 10 15]). 0 reverts nothing, reproducing kcat's own
%     score.
%
% Returns
% -------
% points : struct array
%     one entry per threshold, in the order given: threshold, kcat (the
%     partially-reverted vector), nChanged, meanDev (mean movement in
%     sigma units), rmse.
%
% See also
% --------
% bestParsimonious, nChanged
p = parseGECKOargs(varargin, { ...
    'thresholds', [0 1 2 3 4 6 8 10 15]});
thresholds = p.thresholds;

kcat0     = kcat0(:);
sigma0log = sigma0log(:);
kcat      = kcat(:);

points = repmat(struct('threshold',[],'kcat',[],'nChanged',[],'meanDev',[],'rmse',[]), ...
    numel(thresholds), 1);
for t = 1:numel(thresholds)
    thr = thresholds(t);
    vec = revertBelow(kcat,kcat0,sigma0log,thr);
    points(t).threshold = thr;
    points(t).kcat       = vec;
    points(t).nChanged   = nChanged(vec,kcat0);
    points(t).meanDev    = mean(movementInSigma(vec,kcat0,sigma0log));
    points(t).rmse       = scoreFcn(vec);
end
end

function dev = movementInSigma(kcat,kcat0,sigma0log)
% Per-kcat displacement from the prior, in prior standard deviations.
dev = abs(log(kcat ./ kcat0)) ./ sigma0log;
end

function vec = revertBelow(kcat,kcat0,sigma0log,threshold)
% Put every kcat that moved less than threshold sigma back on its prior.
keep = movementInSigma(kcat,kcat0,sigma0log) >= threshold;
vec = kcat0;
vec(keep) = kcat(keep);
end
