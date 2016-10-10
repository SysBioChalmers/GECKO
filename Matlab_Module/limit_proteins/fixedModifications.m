
function model = fixedModifications(model,options)

%Current values (aerobic Yeast 7.6)
if nargin < 2
    options = [2;3;0.7;1;2;2;0.4266];
end

%Open or close other exchange reactions:
if options(1) == 1  %Open
    model.ub(strcmp(model.rxnNames,'acetate exchange'))              = Inf;
    model.ub(strcmp(model.rxnNames,'formate exchange'))              = Inf;
    model.ub(strcmp(model.rxnNames,'pyruvate exchange'))             = Inf;
    model.ub(strcmp(model.rxnNames,'acetaldehyde exchange'))         = Inf;
    model.ub(strcmp(model.rxnNames,'(R,R)-2,3-butanediol exchange')) = Inf;
else                %Close
    model.ub(strcmp(model.rxnNames,'formate exchange'))              = 1e-2;
    model.ub(strcmp(model.rxnNames,'pyruvate exchange'))             = 5e-2;
%     model.ub(strcmp(model.rxnNames,'acetaldehyde exchange'))         = 0;
%     model.ub(strcmp(model.rxnNames,'(R,R)-2,3-butanediol exchange')) = 0;
end

%Add NGAM reaction:
%           ATP  +  H2O  ->  ADP  +   H+   +  PO4
mets      = {'s_0434','s_0803','s_0394','s_0794','s_1322'};
coefs     = [-1,-1,1,1,1];
[model,~] = addReaction(model,'NGAM',mets,coefs,false,options(3),options(3));

end