%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [ecModel_batch,OptSigma] = getConstrainedModel(ecModel,sigma,Ptot,obj_Val,name)
%
% Function that gets a GEM with kinetic data and returns an enzyme 
% constrained model, either with individual enzyme levels or with the total
% measured protein content.
% 
% The model parameters are automatically curated (querying the BRENDA
% files) for reaching the experimental maximal growth rate on glucose
% minimal media (if it was overconstrained). Then, the average saturation 
% factor for enzymes is fitted for the same growth conditions. Finally the 
% fitted model is simulated on the same conditions and the top ten used 
% enzymes (mass-wise) are saved in a file as an output, this can be used for 
% suggesting further parameters curation (enzyme usages > 10% of the total 
% proteome).
%
% Benjamin Sanchez      2018-08-10
% Ivan Domenzain        2018-09-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel_batch,OptSigma] = getConstrainedModel(ecModel,c_source,sigma,Ptot,obj_Val,modifications,name)
	
    %Get f (estimated mass fraction of enzymes in model)
    [f,~] = measureAbundance(ecModel.enzymes);
    %Get a preliminary enzyme constrained model for performing the Kcats
    %sensitivity analysis
    [ecModel_batch,~,~] = constrainEnzymes(ecModel,Ptot,sigma,f);
	solution            = solveLP(ecModel_batch,1);
    if ~isempty(solution.f)
        %Set the media according to the experimental conditions
        cd ../kcat_sensitivity_analysis
        [ecModel_batch,~] = changeMedia_batch(ecModel_batch,c_source,'Min');
        ObjIndex          = find(ecModel_batch.c);
        % If the model is overconstrained
        if (obj_Val-solution.x(ObjIndex))>0 
            fprintf('\n')
            disp('***************************************************************')
            disp('                The ECmodel is overconstrained                 ')
            %Perform a sensitivity analysis on the objective function with 
            %respect to the individual Kcat coefficients, the algorithm will 
            %iterate replacing the top limiting value according to the maximum 
            %value available in BRENDA for the same EC number until the objective
            %is no longer underpredicted 
            ecModel_batch = modifyKcats(ecModel_batch,obj_Val,modifications,name);
        else
            fprintf('\n')
            disp('***************************************************************')
            disp('              The ECmodel is not overconstrained               ')
        end    
        %The sigma factor is reffited for the specified conditions (constraints in the model)
        fprintf('\n')
        disp('***************************************************************')
        disp('        Fitting the average enzymes saturation factor          ')
        [ecModel_batch,~] = changeMedia_batch(ecModel_batch,c_source,'Min');
        OptSigma          = sigmaFitter(ecModel_batch,Ptot,obj_Val,f);
        enzymePos         = strcmp(ecModel_batch.rxns,'prot_pool_exchange');
        currentEnzymeUB   = ecModel_batch.ub(enzymePos);
        ecModel_batch     = setParam(ecModel_batch,'ub','prot_pool_exchange', ...
                                     currentEnzymeUB*OptSigma/sigma);
        
        %Simulate growth on minimal media and export the top ten used 
        %enzymes to the file "topUsedEnzymes.txt" in the containing folder
        [ecModel_batch,~] = changeMedia_batch(ecModel_batch,c_source,'Min');
        solution          = solveLP(ecModel_batch,1);
        topUsedEnzymes(solution.x,ecModel_batch,'Min_glucose',name);
        cd ../limit_proteins
    else
        disp('ecModel with enzymes pool constraint is not feasible')
    end
end
