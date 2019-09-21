function results = run_ecFSEOF(model,rxnTarget,cSource,alphaLims,Nsteps,file1,file2)
% 
% run_ecFSEOF
%
% 	Function that runs Flux-scanning with Enforced Objective Function
% 	for a specified production target.
%   
%		model     (struct) ecModel with total protein pool constraint
%		rxnTarget (string) Rxn ID for the production target reaction, 
% 				  a exchange reaction is recommended.
%		cSource	  (string) Rxn name for the main carbon source uptake 
%			      reaction
% 		alphaLims (vector) Minimum and maximum biomass yield [gDw/mmol Csource] 
%				  for enforced objective limits
%       Nsteps    (integer) Number of steps for suboptimal objective 
%                 in FSEOF
%		file1 	  (string) opt, file name for an output .txt tab-separated 
% 				  file for the results at the genes level
%		file2 	  (string) opt, file name for an output .txt tab-separated 
% 				  file for the results at the reactions level		
%
% Usage: run_ecFSEOF(ecModel,rxnTarget,cSource,alphaLims,Nsteps)
%
% Last modified.  Ivan Domenzain 2019-09-13
%

if nargin < 7
	file2 = [];
	if nargin < 6
		file1 = [];
	end
end
current       = pwd;
results.model = model;
%Define alpha vector for suboptimal enforced objective values
alphaV        = alphaLims(1):((alphaLims(2)-alphaLims(1))/(Nsteps-1)):alphaLims(2);
%Standardize grRules and rxnGeneMat in model
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules        = grRules;
model.rxnGeneMat     = rxnGeneMat;
%run FSEOF analysis
results.substrate = compare_substrate(model,rxnTarget,cSource,alphaV,1E-4,true);
%Create gene table:
results.geneTable      = cell(length(results.substrate.genes),3);
results.geneTable(:,1) = results.substrate.genes;
results.geneTable(:,2) = results.substrate.geneNames;
results.geneTable(:,3) = num2cell(results.substrate.k_genes);
try
    if ~isempty(file1)
    	T = cell2table(results.geneTable);
    	writetable(T,'Delimiter','\t','QuoteStrings',false) 
    end
    if ~isempty(file2)
    	T = cell2table(results.k_rxns);
    	writetable(T,'Delimiter','\t','QuoteStrings',false) 
    end
catch
    warning('Output files directory not found');
end
end