function [GAMnonPol, GAMpol, GAEC] = splitGAEC(model)
% splitGAM
%
% Splits the growth associated energy cost (GAEC) from the model into:
% 1.    GAMnonPol   the part of GAEC that is an amalgation of non-specified
%                   growth associated maintenance costs
% 2.    GAMpol      the part of GAEC that can be attributed to the cost of
%                   polymerizing macromolecules
% 
% GAEC is taken from the biomass reaction, while the composition of
% macromolecules is used to calculate GAMpol.
%
% Usage: [GAMnonPol, GAMpol, GAEC] = splitGAEC(model)
%

cd ..
parameters = getModelParameters;
bioRxn = parameters.bioRxn;
cost = parameters.pol_cost;
cd limit_proteins

[~,P,C,R,D] = sumBioMass(model);

GAMpol = P*cost(1) + C*cost(2) + R*cost(3) + D*cost(4);

bioRxn = strcmp(model.rxns,bioRxn);
ADP = find(strcmp(model.metNames,'ADP'));
GAMtotal = max(full(model.S(ADP,bioRxn)));

GAMnonPol = GAMtotal - GAMpol;
