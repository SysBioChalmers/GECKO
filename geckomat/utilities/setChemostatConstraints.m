function model = setChemostatConstraints(model,pos,Drate,minProt,flexF,GUR)
% setChemostatConstraints
%
% Function that takes a metabolic model (ecModel, ecModel+Prots or a
% classic GEM) and set chemostat constraints for performing FBA
% simulations.
%
%   model      (struct) MATLAB structure for the metabolic model
%   pos        (vector) Indexes for the main carbon source uptake rate rxn  
%              and the biomass production, respectively.
%   Drate      Dilution rate for the desired chemostat conditions.
%   minProt    TRUE if minimization of total protein usage is desired to be
%              an objective for minimization after a minimal GUR has been 
%              fixed.Default (false)
%   flexF      Relatibve flexibility factor for fixing bounds in the model. 
%              Default (0.001).
%   GUR        Optional. GUR experimental value. It should be provided if
%              it is desired as a fixed condition.
%
% model = setChemostatConstraints(model,pos,Drate,minProt,flexF,GUR)
%
% Last modified.  Ivan Domenzain 2019-09-10

if nargin<6
    GUR = [];
    if nargin<5
        flexF = 0.001;
        if nargin<4
            minProt = false;
        end
    end
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
model.lb(gPos) = (1-flexF)*Drate;
model.c(:)     = 0;
%Set minimization of carbon source uptake as an objective function
if ecFlag
    model.c(cSource)  = -1;
    model.lb(cSource) = 0;
    model.ub(cSource) = 1000;
    if ~isempty(GUR)
        model.lb(cSource) = (1-flexF)*GUR;
    end
else
    model.c(cSource)  = 1;
    model.lb(cSource) = -1000;
    model.ub(cSource) = 0;
    if ~isempty(GUR)
        model.ub(cSource) = (1+flexF)*GUR;
    end
end
%Get GUR minimization solution
solution = solveLP(model);
solution = solution.x;
%If feasible then then fix optimal GUR
if ~isempty(solution)
    if minProt 
        model.lb(cSource) = solution(cSource);
        model.ub(cSource) = (1+flexF)*solution(cSource);
        %Set protein usage as objective to minimize
        model.c(:)       = 0;
        model.c(Protpos) = -1;
    end
else
    disp('Chemostat conditions too stringent')
end
end