function solution = simulateChemostat(model,Drate,pos,minProt)
% setChemostatConstraints
%
% Function that takes a metabolic model (ecModel, ecModel+Prots or a
% classic GEM) and set chemostat constraints for performing FBA
% simulations.
%
%   model      (struct) MATLAB structure for the metabolic model
%   pos        (vector) Indexes for the main carbon source uptake rate rxn  
%              and the biomass procuction one.
%   Drate      Dilution rate for the desired chemostat conditions.
%   minProt    TRUE if minimization of total protein usage is desired to be
%              an objective for minimization after a minimal GUR has been 
%              fixed.
%
% Usage: solution = simulateChemostat(model,Drate,pos,minProt)
%
% Last modified.  Ivan Domenzain 2019-04-15

if nargin<4
  minProt = false;
end
cSource = pos(1);
gPos    = pos(2);
%Identify if model is enzyme-constrained
Protpos = find(contains(model.rxnNames,'prot_pool'));  
if ~isempty(Protpos) 
    ecFlag = true;
else
    ecFlag = false;
end
%Fix dilution rate 
model.lb(gPos) = 0.99*Drate;
model.ub(gPos) = 1.01*Drate;
model.c(:)     = 0;
%Set minimization of carbon source uptake as an objective function
if ecFlag
	model.c(cSource)  = -1;
else
    model.c(cSource)  = 1;
    model.lb(cSource) = -1000;
end
%Run optimization
solution = solveLP(model,1);
solution = solution.x;
if ~isempty(solution)
    if minProt 
        model.lb(cSource) = 0.99999*solution(cSource);
        model.ub(cSource) = 1.00001*solution(cSource);
        %Set protein usage as objective to minimize
        model.c(:)       = 0;
        model.c(Protpos) = -1;
        solution = solveLP(model,1);
        if ~isempty(solution.x)
            solution = solution.x;
        else
            disp('Protein minimization not feasible')
        end
    end
else
    disp('Chemostat conditions too stringent')
end
end