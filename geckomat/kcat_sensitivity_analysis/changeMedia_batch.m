function model = changeMedia_batch(model,c_source,flux)
%changeMedia_batch
%   function that modifies the ecModel and makes it suitable for batch growth
%   simulations on different carbon sources.
%
%   model       (struct) An enzyme constrained model
%   c_source    (string) Rxn name for the main carbon source uptake reaction
%   flux        (doule) Experimental flux value [mmol/gDw h] for the main 
%               carbon source uptake reaction.
%
%
%   Usage: model = changeMedia_batch(model,c_source,flux)
%
% Benjamin J. Sanchez	2018-12-11
% Ivan Domenzain        2019-10-09

% Give the carbon source (c_source) input variable with the following
% format: c_source  = 'D-glucose exchange (reversible)'
if nargin<3
    flux = 1000;
end
%first block any uptake
[rxnIDs,exchange]  = getExchangeRxns(model);
exchange           = exchange(find(contains(rxnIDs,'_REV')));
model.ub(exchange) = 0;
%Allow main carbon source uptake
c_id  = model.rxns(strcmp(model.rxnNames,c_source));
model = setParam(model,'ub',c_id,flux);
%block glucose and oxygen production
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;
%Allow uptake of essential components
model = setParam(model, 'ub', 'r_1654_REV', Inf); % 'ammonium exchange';
model = setParam(model, 'ub', 'r_2100_REV', Inf); % 'water exchange' ;
model = setParam(model, 'ub', 'r_1861_REV', Inf); % 'iron(2+) exchange';
model = setParam(model, 'ub', 'r_1992_REV', Inf); % 'oxygen exchange';
model = setParam(model, 'ub', 'r_2005_REV', Inf); % 'phosphate exchange';
model = setParam(model, 'ub', 'r_2060_REV', Inf); % 'sulphate exchange';
model = setParam(model, 'ub', 'r_1832_REV', Inf); % 'H+ exchange' ;
model = setParam(model, 'ub', 'r_4593_REV', Inf); % 'chloride exchange' ;
model = setParam(model, 'ub', 'r_4595_REV', Inf); % Mn(2+) exchange
model = setParam(model, 'ub', 'r_4596_REV', Inf); % Zn(2+ exchange
model = setParam(model, 'ub', 'r_4597_REV', Inf); % Mg(2+) exchange
model = setParam(model, 'ub', 'r_2049_REV', Inf); % sodium exchange
model = setParam(model, 'ub', 'r_4594_REV', Inf); % Cu(2+) exchange
model = setParam(model, 'ub', 'r_4600_REV', Inf); % Ca(2+) exchange
model = setParam(model, 'ub', 'r_2020_REV', Inf); % potassium exchange
%Block some production fluxes
model = setParam(model, 'ub', 'r_1663', 0); % bicarbonate exchange
model = setParam(model, 'ub', 'r_4062', 0); % lipid backbone exchange
model = setParam(model, 'ub', 'r_4064', 0); % lipid chain exchange
%Allow biomass production 
model = setParam(model, 'ub', 'r_2111', Inf); % growth
end