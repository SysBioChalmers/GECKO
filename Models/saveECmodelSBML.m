%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = saveECmodelSBML(model,toolbox,name,version)
%
% Benjamín J. Sánchez. Last edited: 2018-08-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveECmodelSBML(model,toolbox,name,version)

%Define file path for storage:
if endsWith(name,'_batch')
    file_name = [name(1:strfind(name,'_batch')-1) '/' name];
else
    file_name = [name '/' name];
end

%Model description:
model.description = [name '_' version];

%Format S matrix: avoid long decimals
for i = 1:length(model.mets)
    for j = 1:length(model.rxns)
        if model.S(i,j) ~= 0
            orderMagn    = ceil(log10(abs(model.S(i,j))));
            model.S(i,j) = round(model.S(i,j),6-orderMagn);
        end
    end
end

%Remove model.fields (added by COBRA functions)
if isfield(model,'rules')
    model = rmfield(model,'rules');
end

%Save model as mat:
S.(name) = model;
save([file_name '.mat'], '-struct', 'S')

%Transform model back to COBRA for saving purposes:
if strcmp(toolbox,'COBRA')
    model = ravenCobraWrapper(model);    
    %Remove fields from COBRA model (temporal):
    model = rmfield(model,'metCharges');
    model = rmfield(model,'metChEBIID');
    model = rmfield(model,'metKEGGID');
    model = rmfield(model,'rxnConfidenceScores');
    model = rmfield(model,'rxnECNumbers');
    model = rmfield(model,'rxnKEGGID');
    model = rmfield(model,'rxnReferences');
    model = rmfield(model,'subSystems');
end

%Save model as sbml and text:
writeCbModel(model,'sbml',[file_name '.xml']);
writeCbModel(model,'text',[file_name '.txt']);

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
copyfile([file_name '.xml'],'backup.xml')
fin  = fopen('backup.xml', 'r');
fout = fopen([file_name '.xml'], 'w');
still_reading = true;
while still_reading
  inline = fgets(fin);
  if ~ischar(inline)
      still_reading = false;
  else
      if ~isempty(regexp(inline,'-00[0-9]','once'))
          inline = strrep(inline,'-00','-0');
      elseif ~isempty(regexp(inline,'-01[0-9]','once'))
          inline = strrep(inline,'-01','-1');
      end
      fwrite(fout, inline);
  end
end
fclose('all');
delete('backup.xml');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%