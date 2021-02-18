function results = run_ecFSEOF(model,rxnTarget,cSource,alphaLims,Nsteps,file1,file2)
% 
% run_ecFSEOF
%
% 	Function that runs Flux-scanning with Enforced Objective Function
% 	for a specified production target.
%   
%		model     (struct) ecModel with total protein pool constraint.
%		rxnTarget (string) Rxn ID for the production target reaction, 
% 				  a exchange reaction is recommended.
%		cSource	  (string) Rxn ID for the main carbon source uptake 
%			      reaction (make sure that the correct directionality is
%			      indicated).
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

if nargin < 7
	file2 = [];
	if nargin < 6
		file1 = [];
	end
end
%Check model format, if COBRA model then convert to RAVEN format
if isfield('model','rules') && ~isfield('model','metComps')
    model = ravenCobraWrapper(model);
end
%Define alpha vector for suboptimal enforced objective values
alphaV  = alphaLims(1):((alphaLims(2)-alphaLims(1))/(Nsteps-1)):alphaLims(2);
%Standardize grRules and rxnGeneMat in model
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules        = grRules;
model.rxnGeneMat     = rxnGeneMat;
%run FSEOF analysis
results = ecFlux_scanning(model,rxnTarget,cSource,alphaV,1E-4,true);
%Create gene table:
results.geneTable      = cell(length(results.genes),3);
results.geneTable(:,1) = results.genes;
results.geneTable(:,2) = results.geneNames;
results.geneTable(:,3) = num2cell(results.k_genes);
%Create rxns table (exclude enzyme usage reactions):
toKeep                 = find(~startsWith(results.rxns(:,1),'draw_prot_'));
results.k_rxns         = results.k_rxns(toKeep);
results.k_matrix       = results.k_matrix(toKeep,:);
results.v_matrix       = results.v_matrix(toKeep,:);
results.rxnsTable      = cell(length(results.k_rxns),5);
results.rxnsTable(:,1) = results.rxns(toKeep,1);
results.rxnsTable(:,2) = results.rxns(toKeep,2);
results.rxnsTable(:,3) = num2cell(results.k_rxns);
results.rxnsTable(:,4) = results.rxns(toKeep,3);
results.rxnsTable(:,5) = results.rxns(toKeep,4);
try
    if ~isempty(file1)
        varNames = {'gene_IDs' 'gene_names' 'K_score'}; 
    	T        = cell2table(results.geneTable,'VariableNames',varNames);
    	writetable(T,file1,'Delimiter','\t','QuoteStrings',false) 
    end
    if ~isempty(file2)
        varNames = {'rxn_IDs' 'rxnNames' 'K_score' 'grRules' 'rxn_formula'};
    	T        = cell2table(results.rxnsTable,'VariableNames',varNames);
        writetable(T,file2,'Delimiter','\t','QuoteStrings',false) 
    end
catch
    warning('Output files directory not found');
end
%Remove redundant output fields
results = rmfield(results,'k_rxns');
results = rmfield(results,'rxns');
results = rmfield(results,'genes');
results = rmfield(results,'geneNames');
results = rmfield(results,'k_genes');
end