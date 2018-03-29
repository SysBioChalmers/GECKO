%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [ecModel_batch,OptSigma] = getConstrainedModel(ecModel,sigma,Ptot,gR_exp,name)
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
% suggesting further parameter curation targets 
% (enzyme usages > 10% of the total proteome).
%
% Ivan Domenzain    Last edited. 2018-03-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel_batch,OptSigma] = getConstrainedModel(ecModel,sigma,Ptot,gR_exp,modifications,name)
	current = pwd;
	%Get a preliminary enzyme constrained model for performing the Kcats 
	%sensitivity analysis
	ecModel_batch = constrainEnzymes(ecModel,Ptot,sigma);
    ObjIndex      = find(ecModel_batch.c);
	solution      = solveLP(ecModel_batch,1);
    if ~isempty(solution.f)
        %Set the media according to the media of the experimental measurement
        cd ../Kcat_sensitivity_analysis
        c_source          = 'D-glucose exchange (reversible)';
        %[ecModel_batch,~] = changeMedia_batch(ecModel_batch,c_source,'Min');
        carbonUptake = 1000;%5.8;%[mmol/gDwh]
        OxyUptake    = 1000;%13; %[mmol/gDwh]
        [ecModel_batch,~] = fixUptakeRates(ecModel_batch,c_source,'Min',[],carbonUptake,OxyUptake);
        solution          = solveLP(ecModel_batch,1);
        % If the model is overconstrained
        if ((gR_exp-solution.x(ObjIndex))/gR_exp)>0.05 
            disp('***************************************************************')
            disp('                The ECmodel is overconstrained                 ')
            %Perform a sensitivity analysis on the individual Kcat coefficients,
            %the algorithm will iterate changing the top growth limiting value 
            %according to the maximum value available in BRENDA for the same 
            %EC number until the growth rate is no longer underpredicted 
            ecModel = modifyKcats(ecModel,ecModel_batch,gR_exp,modifications,name);
        else
            disp('***************************************************************')
            disp('              The ECmodel is not overconstrained               ')
        end
        %The sigma factor is reffited for minimal glucose media
        disp('***************************************************************')
        disp('        Fitting the average enzymes saturation factor          ')
        OptSigma = sigmaFitter(ecModel,Ptot,gR_exp);
        %Avoid extreme values
        if OptSigma == 0.4 || OptSigma == 0.6
            OptSigma = 0.5;
        end
        %The ecModel with new modified Kcat values is constrained with the 
        %optimal sigma value found
        cd ../limit_proteins
        ecModel_batch = constrainEnzymes(ecModel,Ptot,OptSigma);
        solution      = solveLP(ecModel_batch,1);
        % If the model is overconstrained
        if ((gR_exp-solution.x(ObjIndex))/gR_exp)>0.05
            disp('***************************************************************')
            disp('             Fitting the ATP growth requeriments               ')
            biomassRxn    = ecModel_batch.rxnNames(ObjIndex-1);
            cd ../Kcat_sensitivity_analysis
            ecModel_batch = fitGAM(ecModel_batch,biomassRxn,gR_exp,Ptot);
        end            
        %Simulate growth on minimal glucose media and export the top ten used 
        %enzymes to the file "topUsedEnzymes.txt" in the containing folder
        cd (current)
        cd ../Kcat_sensitivity_analysis
        c_source          = 'D-glucose exchange (reversible)';
        [ecModel_batch,~] = changeMedia_batch(ecModel_batch,c_source,'Min');
        solution          = solveLP(ecModel_batch,1);
        topUsedEnzymes(solution.x,ecModel_batch,'Min_glucose',name);
    else
        disp('ecModel with enzymes pool constraint is not feasible')
    end
end