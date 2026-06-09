function [model, sigma] = sigmaFitter(model, growthRate, Ptot, f, makePlot, modelAdapter)
% sigmaFitter  Fit the average enzyme saturation factor of an ecModel.
%
% Fits the average enzyme saturation factor in an ecModel according to a
% provided experimentally measured value for the objective function (i.e.
% growth rate at specified conditions).
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% growthRate : double, optional
%     growth rate that should be reached (default: read from the model
%     adapter).
% Ptot : double, optional
%     total cellular protein content in g/gDCW (default: read from the model
%     adapter; if not specified there, 0.5 g/gDCW is assumed).
% f : double, optional
%     estimated fraction of enzymes in the model (default: read from the
%     model adapter; if not specified there, 0.5 is assumed).
% makePlot : logical, optional
%     whether a plot should be made (default true).
% modelAdapter : ModelAdapter, optional
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
%     [model, sigma] = sigmaFitter(model, growthRate, Ptot, f, makePlot, modelAdapter);

if nargin < 6 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if nargin<5 || isempty(makePlot)
    makePlot = true;
end
if nargin<4 || isempty(f)
    f = modelAdapter.getParameters().f;
end
if nargin<3 || isempty(Ptot)
    Ptot = modelAdapter.getParameters().Ptot;
end
if nargin<2 || isempty(growthRate)
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
if makePlot
    figure
    plot(sigParam,errors,'LineWidth',5)
    title('Sigma fitting')
    xlabel('Average enzyme saturation [-]')
    ylabel('Absolute relative error [%]')
end
end
