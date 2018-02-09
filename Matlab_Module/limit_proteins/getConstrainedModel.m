%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ecModel_batch = getConstrainedModel(ecModel);
%
% Function that gets a model with kinetic data and returns an enzyme 
% constrained model, either with individual enzyme levels or with the total
% measured protein content.
%
% Ivan Domenzain    Last edited. 2018-02-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ecModel_batch = getConstrainedModel(ecModel,sigma,Ptot,gR_exp)
	
	%Get a preliminary enzyme constrained model for performing the Kcats 
	%sensitivity analysis
	ecModel_batch = constrainEnzymes(ecModel,Ptot,sigma);
	
	%Set the media according to the media of the experimental measurement
	cd ../Kcat_sensitivity_analysis
	c_source          = 'D-glucose exchange (reversible)';
	[ecModel_batch,~] = changeMedia_batch(ecModel_batch,c_source,'Min');

	%Perform a sensitivity analysis on the individual Kcat coefficients,
	%the algorithm will iterate changing the top growth limiting value 
	%according to the maximum value available in BRENDA for the same 
	%EC number until the growth rate is no longer underpredicted 
	ecModel           = modifyKcats(ecModel,ecModel_batch,0.41);
	
	%The sigma factor is reffited for minimal glucose media 
	OptSigma          = sigmaFitter(ecModel,Ptot,gR_exp);
	
	%The ecModel with new modified Kcat values is constrained with the 
	%optimal sigma value found
	ecModel_batch     = constrainEnzymes(ecModel,Ptot,OptSigma);
	
	%Simulate growth on minimal glucose media and export the top ten used 
	%enzymes to the file "topUsedEnzymes.txt" in the containing folder
	cd ../Kcat_sensitivity_analysis
	solution = solveLP(ecModel_batch);
	topUsedEnzymes(solution.x,ecModel_batch,'Min_glucose');
end