%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = saveECmodelSBML(model,name,isBatch)
%
% Benjamín J. Sánchez. Last edited: 2018-03-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveECmodelSBML(model,name,isBatch)

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

%Format gene rule field:
for i = 1:length(model.rxns)
    if ~isempty(model.rules{i})
        pos = find(strcmp(model.genes,model.rules{i}));
        if ~isempty(pos)
            model.rules{i} = ['x(' num2str(pos) ')'];
        end
    end
end
model = rmfield(model,'grRules');

%Format metFormulas:
model.metFormulas = strrep(model.metFormulas,'(','');
model.metFormulas = strrep(model.metFormulas,')n','');
model.metFormulas = strrep(model.metFormulas,')','');

%Batch case: modify name
folder = name;
if isBatch
    name = [name '_batch'];
end

%Save model:
writeCbModel(model,'sbml',[folder '/' name '.xml']);
writeCbModel(model,'text',[folder '/' name '.txt']);

%Remove lines of the sort "<fbc:geneProductAssociation/>" from xml file, and
%convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
copyfile([folder '/' name '.xml'],'backup.xml')
fin  = fopen('backup.xml', 'r');
fout = fopen([folder '/' name '.xml'], 'w');
still_reading = true;
while still_reading
  inline = fgets(fin);
  if ~ischar(inline)
      still_reading = false;
  elseif ~contains(inline,'<fbc:geneProductAssociation/>')
      inline = strrep(inline,'-00','-0');
      fwrite(fout, inline);
  end
end
fclose('all');
delete('backup.xml');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%