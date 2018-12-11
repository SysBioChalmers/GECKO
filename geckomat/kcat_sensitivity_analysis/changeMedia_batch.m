%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeMedia(model,media,flux)
%
% function that modifies the ecModel and makes it suitable for batch growth
% simulations on different carbon sources.
%
% INPUT:
%   - model:  An enzyme constrained model
%   - meadia: Media type ('YEP' for complex, 'MAA' minimal with Aminoacids,
%                          'Min' for minimal media)
%   - flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh
%
% OUTPUT:
%   - model: The ECmodel with
%
% Ivan Domenzain        2018-09-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeMedia_batch(model,c_source)

c_id  = model.rxns(strcmp(model.rxnNames,c_source));
model = setParam(model,'ub',c_id,Inf);
