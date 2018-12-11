%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeMedia_batch(model,c_source)
%
% Function that modifies the ecModel and makes it suitable for batch growth
% simulations on the carbon source of choice. You can add more changes in
% this function if you need to do so for your model.
%
% INPUT:
%   model       An enzyme constrained model.
%	c_source	The name of the exchange reaction that supplies the model
%				with carbon.
%
% OUTPUT:
%   model       The enzyme constrained model with modified boundaries.
%
% Ivan Domenzain        2018-09-27
% Benjamin J. Sanchez	2018-12-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeMedia_batch(model,c_source)

c_id  = model.rxns(strcmp(model.rxnNames,c_source));
model = setParam(model,'ub',c_id,Inf);
