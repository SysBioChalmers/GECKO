%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OptSigma = sigmaFitter(model_batch,Ptot,expVal,f)
% 
% Function that fits the average enzyme saturation factor in an ecModel
% according to a provided experimentally measured value for the objective
% function (i.e. growth rate at specified conditions)
%
% INPUTS:
%       model_batch     An EC batch model with an initial sigma factor
%                       assigned
%       Ptot            Total protein amount in the model (Experimental)
%                       [g/gDw]
%       expVal          Experimentally measured value for the objective function
%       f               Estimated mass fraction of enzymes in model [g/g]
%
% OUTPUTS:
%       optSigma    The optimal sigma value obtained
%
% Benjamin Sanchez      2018-08-10
% Ivan Domenzain        2018-09-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptSigma = sigmaFitter(model_batch,Ptot,expVal,f)

objValues   = [];
errors    = [];
sigParam  = [];
poolIndex = find(strcmpi(model_batch.rxnNames,'prot_pool_exchange'));
objPos    = find(model_batch.c);
%Relax bounds for the objective function
model_batch.lb(objPos) = 0;
model_batch.ub(objPos) = 1000;

for i=1:100
    %Constrains the ecModel with the i-th sigma factor
    sigma = i/100;
    model_batch.ub(poolIndex) = sigma*Ptot*f; 
    solution  = solveLP(model_batch,1);
    if isempty(solution.x)
        solution.x=zeros(length(model_batch.rxns),1);
    end
    objValues = [objValues; solution.x(objPos)];
    error     = abs(((expVal-solution.x(objPos))/expVal)*100);
    errors    = [errors; error];
    error     = num2str(((expVal-solution.x(objPos))/expVal)*100);
    sigParam  = [sigParam; sigma];
    disp(['Fitting sigma factor: ' num2str(sigma) '   error: ' error '%'])
end
[~, minIndx] = min(errors);
OptSigma     = sigParam(minIndx);
figure
plot(sigParam,errors,'LineWidth',5)
title('Sigma fitting for growth on glucose minimal media')
xlabel('Average enzyme saturation [-]')
ylabel('Absolute relative error [%]')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
