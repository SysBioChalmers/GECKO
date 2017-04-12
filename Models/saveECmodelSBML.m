%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = saveECmodelSBML(model,name)
%
% Benjamín J. Sánchez. Last edited: 2017-04-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = saveECmodelSBML(model,name)

%Introduce compartments to both metabolite ID and name:
comps     = model.comps;
compNames = model.compNames;
for i = 1:length(model.mets)
    comp_ID           = comps{model.metComps(i)};
    comp_name         = compNames{model.metComps(i)};
    model.mets{i}     = [model.mets{i} '[' comp_ID ']'];
    model.metNames{i} = [model.metNames{i} ' [' comp_name ']'];
end

%Remove unused fields:
model = rmfield(model,'metComps');
model = rmfield(model,'comps');
model = rmfield(model,'compNames');

%Fix gene names;
model.genes   = regexprep(model.genes,'-','_');
model.grRules = regexprep(model.grRules,'-','_');

%Save model:
model.id = 'ecYeast7';
writeCbModel(model,'sbml',[name '.xml'],comps,compNames,3,1);
writeCbModel(model,'text',[name '.txt'],comps,compNames,3,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%