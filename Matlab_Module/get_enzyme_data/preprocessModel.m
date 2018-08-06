%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = preprocessModel(model)
% Deletes blocked rxns + corrects model.rev
%
% INPUT:    A model structure
% OUTPUT:   The corrected model
%
% Benjamín J. Sánchez. Last edited: 2018-08-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = preprocessModel(model)

%Delete blocked rxns (LB = UB = 0):
to_remove = boolean((model.lb == 0).*(model.ub == 0));
model     = removeRxns(model,model.rxns(to_remove));

%Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%