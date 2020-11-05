function [model,NGAM] = fitNGAM(model,NGAMrxn,exp_data,bounds,ECflag,Drate)
% fitNGAM
%
%  Function that takes an ecModel for a given organism and fits the LB for 
%  the non-growth associated maintenance pseudoreaction (NGAM) reaction 
%  through an error minimization with respect to the provided experimental 
%  data.
% 
%     model    GEM or ecGEM matlab structure.
%     NGAMrxn  rxnID for the NGAM reaction.
%     exp_data (vector) Experimental values for GUR O2_EX and CO2_Ex [mmol/gDw h]
%     bounds   (vector) Prior bounds for the NGAM value [min max].
%     ECflag   (logical) True if model is in irreversible format
%     Drate    (optional) Dilution rate if chemostat conditions are desired
%              to be set prior to the fitting process.
%
% Created.       Ivan Domenzain 2019-02-21
% Last modified. Ivan Domenzain 2019-04-24

current = pwd;
if nargin<5
    ECflag = true;
end
NGAMIndex           = find(strcmpi(model.rxns,NGAMrxn));
model.lb(NGAMIndex) = bounds(1);
model.ub(NGAMIndex) = +1000;
%Relevant positions:
if ECflag 
    cd ../..
    parameters = getModelParameters;
    c_source   = parameters.exch_names{2};
    oxy_exch   = parameters.exch_names{3};
    CO2_exch   = parameters.exch_names{4};
    pos(1)     = find(strcmp(model.rxnNames,c_source));
    pos(2)     = find(strcmp(model.rxnNames,CO2_exch));
    pos(3)     = find(strcmp(model.rxnNames,oxy_exch));
    minProt    = true;
else
    pos(1) = find(strcmp(model.rxnNames,'D-glucose exchange'));
    pos(2) = find(strcmp(model.rxnNames,'carbon dioxide exchange'));
    pos(3) = find(strcmp(model.rxnNames,'oxygen exchange'));
    minProt = false;
end
%Set chemostat conditions (optional)
if nargin>5
	cd ..
    indexes = [pos(1) find(model.c)];
	model   = setChemostatConstraints(model,indexes,Drate,minProt);
end
cd (current)
%Solve
steps = 100;
for i=1:steps
    delta   = (bounds(2)-bounds(1))*(i-1)/steps;
    NGAM(i) = bounds(1) + delta;
    model.lb(NGAMIndex) = NGAM(i);
    sol   = solveLP(model,1);
    %Store relevant variables:
    if ~isempty(sol.x)
        mod_data = sol.x(pos)'; 
    else
        mod_data = zeros(1,length(pos));
    end
    Residues   = (abs(mod_data) - exp_data)./exp_data;
    fitting(i) = sqrt(sum(sum(Residues.^2)));
end
%Choose best fit:
[minErr,best] = min(fitting);

if best == length(NGAM)
    warning('NGAM found is sub-optimal: please expand GAM search bounds.')
    NGAM = bounds(1);
else
    NGAM = NGAM(best);
end

model.lb(NGAMIndex) = NGAM;
disp(['Fitted NGAM = ' num2str(NGAM) ' [mmol ATP/gDw h] -> Error = ' num2str(minErr)])
end