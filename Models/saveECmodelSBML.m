%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = saveECmodelSBML(model,name)
%
% Benjamín J. Sánchez. Last edited: 2017-10-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveECmodelSBML(model,name)

%Introduce compartments to both metabolite ID and name:
comps     = model.comps;
compNames = model.compNames;
for i = 1:length(model.mets)
    comp_ID           = comps{model.metComps(i)};
    comp_name         = compNames{model.metComps(i)};
    model.mets{i}     = [model.mets{i} '[' comp_ID ']'];
    model.metNames{i} = [model.metNames{i} ' [' comp_name ']'];
end
model = rmfield(model,'metComps');

%Format gene field:
model.genes = strrep(model.genes,'-','_');
for i = 1:length(model.rxns)
    if ~isempty(model.rules{i})
        model.rules{i} = strrep(model.rules{i},'-','_');
    end
end
model = rmfield(model,'grRules');

%Format metFormulas:
model.metFormulas = strrep(model.metFormulas,'(','');
model.metFormulas = strrep(model.metFormulas,')n','');
model.metFormulas = strrep(model.metFormulas,')','');

%Save model:
model.id = name;
writeCbModel(model,'sbml',[name '.xml']);
writeCbModel(model,'text',[name '.txt']);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%