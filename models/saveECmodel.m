%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = saveECmodel(model,toolbox,name,version)
%
% Benjamin J. Sanchez. Last edited: 2018-10-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = saveECmodel(model,toolbox,name,version)

%Define file path for storage:
struct_name = 'ecModel';
if endsWith(name,'_batch')
    struct_name = [struct_name '_batch'];
    root_name   = name(1:strfind(name,'_batch')-1);
else
    root_name = name;
end
file_name = [root_name '/' name];

%Model description:
model.description = [struct_name ' of ' lower(root_name(3)) root_name(4:end)];
model.id          = [name '_v' version];

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

if strcmp(toolbox,'COBRA')
    %Transform model back to COBRA for saving purposes:
    model_cobra = ravenCobraWrapper(model);
    %Remove fields from COBRA model (temporal):
    model_cobra = rmfield(model_cobra,'metCharges');
    model_cobra = rmfield(model_cobra,'metChEBIID');
    model_cobra = rmfield(model_cobra,'metKEGGID');
    model_cobra = rmfield(model_cobra,'metSBOTerms');
    model_cobra = rmfield(model_cobra,'rxnConfidenceScores');
    model_cobra = rmfield(model_cobra,'rxnECNumbers');
    model_cobra = rmfield(model_cobra,'rxnKEGGID');
    model_cobra = rmfield(model_cobra,'rxnReferences');
    model_cobra = rmfield(model_cobra,'subSystems');
    model_cobra = rmfield(model_cobra,'rxnSBOTerms');
    %Save model as sbml and text:
    writeCbModel(model_cobra,'sbml',[file_name '.xml']);
    writeCbModel(model_cobra,'text',[file_name '.txt']);
else
    exportForGit(model,name,root_name,{'xml','yml','txt','mat'});
end

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
