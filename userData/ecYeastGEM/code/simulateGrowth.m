function flux = simulateGrowth(model,target,C_source,objCoeff,alpha,tol)
%simulateGrowth
%
%   Function that performs a series of LP optimizations on an ecModel,
%   by first maximizing biomass, then fixing a suboptimal value and 
%   proceeding to maximize a given production target reaction, last 
%   a protein pool minimization is performed subject to the optimal 
%   production level.
%   
%       model     (struct) ecModel with total protein pool constraint.
%                 the model should come with growth pseudoreaction as 
%                 an objective to maximize.
%       target    (string) Rxn ID for the objective reaction.
%       cSource   (string) Rxn ID for the main carbon source uptake 
%                 reaction
%       objCoeff  (integer) Coefficient for the target reaction in the
%                 objective function
%       alpha     (dobule) scalling factor for desired suboptimal growth
%       tol       (double) numerical tolerance for fixing bounds
%
% Usage: flux = simulateGrowth(model,target,C_source,alpha,tol)
%

if nargin<6
    tol = 1E-6;
    if nargin<5
        alpha = 1;
        if nargin<4
            objCoeff = 1;
        end
    end
end
%Fix a unit main carbon source uptake
cSource_indx           = find(strcmpi(model.rxns,C_source));
model.lb(cSource_indx) = (1-tol);
model.ub(cSource_indx) = (1+tol);
%Position of target rxn:
posP      = strcmpi(model.rxns,target);
growthPos = find(model.c);
%Max growth:
sol = solveLP(model);
%Fix growth suboptimal and then max product:
model.lb(growthPos) = sol.x(growthPos)*(1-tol)*alpha;
flux                = optModel(model,posP,objCoeff,tol,sol.x);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = optModel(model,pos,c,tol,base_sol)
%Find reversible reaction for the production target, if available 
if any(pos)
    rxn_code = model.rxns{pos};
    if endsWith(rxn_code,'_REV')
        rev_pos = strcmp(model.rxns,rxn_code(1:end-4));
    else
        rev_pos = strcmp(model.rxns,[rxn_code '_REV']);
    end
    %Block reversible production target reaction
    model.lb(rev_pos) = 0;
    model.ub(rev_pos) = 0;
end
%Optimize for a given objective
model.c      = zeros(size(model.rxns));
model.c(pos) = c;
sol          = solveLP(model);
%If not successful then relax reversible rxn bounds with values
%from the basal distribution
if isempty(sol.x) && ~isempty(base_sol)
    model.lb(rev_pos) = base_sol(rev_pos);
    model.ub(rev_pos) = base_sol(rev_pos);
    sol               = solveLP(model);
end
%Now fix optimal value for objective and minimize unmeasured proteins
%usage
if ~isempty(sol.x)
    model.lb(pos)    = (1-tol)*sol.x(pos);
    protPos          = contains(model.rxnNames,'prot_pool_');
    model.c(:)       = 0;
    model.c(protPos) = -1;
    sol              = solveLP(model,1);
    sol              = sol.x;
else
    sol = zeros(length(model.rxns),1);
end
end