%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = saveECmodelSBML(model,name)
%
% Benjamín J. Sánchez. Last edited: 2017-10-25
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

if nargin == 2
    %Reformat gene structure:
    model.genes = strrep(model.genes2,'-','_');
    model.rules = strrep(model.grRules,'-','_');
	model       = rmfield(model,'genes2');
    model       = rmfield(model,'grRules');
    
    %Reformat metFormulas:
    model.metFormulas = strrep(model.metFormulas,'(','');
    model.metFormulas = strrep(model.metFormulas,')n','');
    model.metFormulas = strrep(model.metFormulas,')','');
    
    %Save model:
    model.id = 'ecYeast7';
    writeCbModel(model,'sbml',[name '.xml'],comps,compNames,3,1);
    writeCbModel(model,'text',[name '.txt'],comps,compNames,3,1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%