function [model, sigma] = sigmaFitter(model, varargin)
% sigmaFitter  Fit the average enzyme saturation factor of an ecModel.
%
% Scans candidate sigma-factors from 0.01 to 1 in steps of 0.01 and picks the
% one whose simulated growth rate best matches the provided experimental
% growth rate, then applies it to the model's protein pool constraint.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% growthRate : double
%     growth rate that should be reached (default: read from the model
%     adapter).
% Ptot : double
%     total cellular protein content in g/gDCW (default: read from the model
%     adapter; if not specified there, 0.5 g/gDCW is assumed).
% f : double
%     estimated fraction of enzymes in the model (default: read from the
%     model adapter; if not specified there, 0.5 is assumed).
% makePlot : logical
%     whether a plot should be made (default true).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% model : struct
%     ecModel with protein pool exchange upper bound adapted to the optimal
%     sigma-factor.
% sigma : double
%     optimal sigma-factor.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, sigma] = sigmaFitter(model);
%     [model, sigma] = sigmaFitter(model, 'makePlot', false);
%     [model, sigma] = sigmaFitter(model, growthRate, Ptot, f, makePlot, modelAdapter);

p = parseGECKOargs(varargin, { ...
    'growthRate',   []; ...
    'Ptot',         []; ...
    'f',            []; ...
    'makePlot',     true; ...
    'modelAdapter', []});
growthRate   = p.growthRate;
Ptot         = p.Ptot;
f            = p.f;
makePlot     = p.makePlot;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if isempty(makePlot)
    makePlot = true;
end
if isempty(f)
    f = modelAdapter.getParameters().f;
end
if isempty(Ptot)
    Ptot = modelAdapter.getParameters().Ptot;
end
if isempty(growthRate)
    growthRate = modelAdapter.getParameters().gR_exp;
end

objValues = [];
errors    = [];
sigParam  = [];
objPos    = find(model.c);
%Relax bounds for the objective function
model.lb(objPos) = 0;
model.ub(objPos) = 1000;
hsSol=[];
for i=1:100
    %Constrains the ecModel with the i-th sigma factor
    sigma = i/100;
    model = setProtPoolSize(model, Ptot, f, sigma, modelAdapter);
    [solution, hsSol]  = solveLP(model,0,[],hsSol);
    if isempty(solution.x)
        solution.x=zeros(length(model.rxns),1);
    end
    objValues = [objValues; solution.x(objPos)];
    error     = abs(((growthRate-solution.x(objPos))/growthRate)*100);
    errors    = [errors; error];
    error     = num2str(((growthRate-solution.x(objPos))/growthRate)*100);
    sigParam  = [sigParam; sigma];
end
[~, minIndx] = min(errors);
sigma     = sigParam(minIndx);
% The loop above leaves the protein pool sized for the *last* sigma tried
% (i=100, i.e. sigma=1), not the best-fitting one found by the search --
% re-apply it here so the returned model actually matches the returned
% sigma, per this function's own documented contract ("ecModel with
% protein pool exchange upper bound adapted to the optimal sigma-factor").
model = setProtPoolSize(model, Ptot, f, sigma, modelAdapter);
if makePlot
    figure
    plot(sigParam,errors,'LineWidth',5)
    title('Sigma fitting')
    xlabel('Average enzyme saturation [-]')
    ylabel('Absolute relative error [%]')
end
end
