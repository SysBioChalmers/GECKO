
function model = fixedModifications(model,chemostat)

%Current values (aerobic Yeast 7.6)
NGAM = 0.7;

%Add NGAM reaction:
%           ATP  +  H2O  ->  ADP  +   H+   +  PO4
mets      = {'s_0434','s_0803','s_0394','s_0794','s_1322'};
coefs     = [-1,-1,1,1,1];
[model,~] = addReaction(model,'NGAM',mets,coefs,false,NGAM,NGAM);

if nargin == 2    
    %Limit measured excretions:
    model.ub(strcmp(model.rxnNames,'acetate exchange'))  = chemostat(1);
    model.ub(strcmp(model.rxnNames,'pyruvate exchange')) = chemostat(2);
    
    %Limit unmeasured excretions:
    model.ub(strcmp(model.rxnNames,'(R,R)-2,3-butanediol exchange')) = 1e-5;
    model.ub(strcmp(model.rxnNames,'acetaldehyde exchange'))         = 1e-5;
    model.ub(strcmp(model.rxnNames,'glycine exchange'))              = 1e-5;
end

end