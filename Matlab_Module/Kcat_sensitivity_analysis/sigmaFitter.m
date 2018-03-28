%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function OptSigma = sigmaFitter(model,Ptot,gR_exp)
% 
% Function that recieves an EC model then finds the top limiting Kcat values
% and tries to relax them according to the available data in BRENDA (Kcats
% and SA*Mw entries).
%
% INPUTS:
%       model       An EC batch model with an initial sigma factor assigned
%       model_data  The model_data structure saved by the GECKO
%                   enhanceGEM.m script
%       Ptot        Total protein amount in the model (Experimental)
%                   [g/gDw]
%       gR_exp      Experimental growth rate on glucose minimal media 
%                   [g/gDw hr]
%
% OUTPUTS:
%       optSigma    The optimal sigma value obtained
%
% Ivan Domenzain.   Last edited 2017-11-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OptSigma = sigmaFitter(model,Ptot,gR_exp)
    gRate_sim = [];error = []; sigParam = [];
    c_source  = 'D-glucose exchange (reversible)';
    [model,~] = changeMedia_batch(model,c_source,'Min');
    cd ../limit_proteins
    for i=40:60
        % Constrains the ecModel with the current sigma factor
        sigma        = i/100;
        %model_batch  = changeCultureMedia(model);
        model_batch  = constrainEnzymes(model,Ptot,sigma);
        % Change to minimal glucose media
        gR_pos        = find(strcmpi(model_batch.rxnNames,'growth'));
        model_batch.c = zeros(size(model_batch.c));
        model_batch.c(gR_pos) = 1;
        solution               = solveLP(model_batch);
        model_batch.lb(gR_pos) = 0.999*solution.x(gR_pos);
        model_batch.ub(gR_pos) = solution.x(gR_pos);
        solution               = solveLP(model_batch,1);
        gRate_sim              = [gRate_sim; solution.x(gR_pos)];
        error                  = [error; abs((gR_exp-solution.x(gR_pos))/gR_exp)*100];
        sigParam               = [sigParam; sigma]; 
        disp(sigma)
        disp(((gR_exp-solution.x(gR_pos))/gR_exp)*100)
    end
    [minError, minIndx] = min(error);
    OptSigma            = sigParam(minIndx);
    figure
    plot(sigParam,error,'LineWidth',5)
    title('Sigma fitting for growth on glucose minimal media')
    xlabel('Average enzyme saturation [-]')
    ylabel('Absolute relative error [%]')
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

