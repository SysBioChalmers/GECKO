%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function model = sensitivityTuning(model, desiredGrowthRate, modelAdapter)
%
% Function that relaxes the most limiting kcats until a certain growth rate 
% is reached.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = sensitivityTuning(model, desiredGrowthRate, modelAdapter)

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

params = modelAdapter.params;


lastGrowth = 0;

m = model;
m.c = double(strcmp(m.rxns, params.bioRxn));%biomass_human set objective function to growth

%To avoid numerical issues, make sure no kcat is below 0.1
%This is not desirable to do, but unfortunately necessary - otherwise the 
%solver says the solution is infeasible
m.ec.kcat(m.ec.kcat < 0.1) = 0.1;
m = applyKcatConstraints(m);

if ~m.ec.geckoLight
    error('not impl')
else
    origRxns = extractAfter(m.ec.rxns,4);
    iteration = 1;
    while lastGrowth < desiredGrowthRate
        protPoolStoich = m.S(strcmp(m.mets, 'prot_pool'),:).';
        res = solveLP(m,0); %skip parsimonius, only takes time
        lastGrowth = -res.f;
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]) 
        iteration = iteration + 1;
        %find the highest protein usage flux
        [~,sel] = min(res.x .* protPoolStoich); %max consumption
        rxn = m.rxns(sel.');
        targetSubRxns = strcmp(origRxns, rxn);
        m.ec.kcat(targetSubRxns) = m.ec.kcat(targetSubRxns) .* 10;
        m = applyKcatConstraints(m,targetSubRxns);
    end
end

model = m;

end
