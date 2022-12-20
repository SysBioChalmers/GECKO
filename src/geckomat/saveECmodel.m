function model = saveECmodel(model,toolbox,name,version,path,condition)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%       model: an ecModel to save
%
%       toolbox: toolbox version to save. 'RAVEN' or 'COBRA'
%
%       name: the file name
%
%       version: model version
%
% Usage:
%
%       model = saveECmodel(model,toolbox,path,name,version)
%
% Benjamin J. Sanchez. Last edited: 2018-10-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['Saving ' name '_' version ':\n'])

% if a condition is defined then save the model as condition
if nargin == 6
    path = [path condition '/'];
    name = condition;
end

%Model description:
model.description = ['Enzyme-constrained model of ' version];
model.id          = [name];

%Format S matrix: avoid long decimals
for i = 1:length(model.mets)
    for j = 1:length(model.rxns)
        if model.S(i,j) ~= 0
            orderMagn    = ceil(log10(abs(model.S(i,j))));
            model.S(i,j) = round(model.S(i,j),6-orderMagn);
        end
    end
end

%For functional models, save upper bounds as +1000:
model.ub(isinf(model.ub)) = 1000;

%Remove model.rules (added by COBRA functions)
model = takeOutField(model,'rules');

if strcmp(toolbox,'COBRA')
    %Transform model back to COBRA for saving purposes:
    model_cobra = ravenCobraWrapper(model);
    %Remove fields from COBRA model (temporal):
    model_cobra = takeOutField(model_cobra,'metCharges');
    model_cobra = takeOutField(model_cobra,'metChEBIID');
    model_cobra = takeOutField(model_cobra,'metKEGGID');
    model_cobra = takeOutField(model_cobra,'metNotes');
    model_cobra = takeOutField(model_cobra,'metSBOTerms');
    model_cobra = takeOutField(model_cobra,'rxnConfidenceScores');
    model_cobra = takeOutField(model_cobra,'rxnECNumbers');
    model_cobra = takeOutField(model_cobra,'rxnKEGGID');
    model_cobra = takeOutField(model_cobra,'rxnReferences');
    model_cobra.subSystems = cell(size(model_cobra.rxns));
    model_cobra = takeOutField(model_cobra,'rxnSBOTerms');
    %Save model as sbml and text:
    writeCbModel(model_cobra,'mat',[path name '.mat']);
    writeCbModel(model_cobra,'sbml',[path name '.xml']);
    writeCbModel(model_cobra,'text',[path name '.txt']);
else
    exportForGit(model,name,path,{'xml','yml','txt','mat'},false,false);
end

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
copyfile([path name '.xml'],'backup.xml')
fin  = fopen('backup.xml', 'r');
fout = fopen([path name '.xml'], 'w');
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

function model = takeOutField(model,field)

if isfield(model,field)
    model = rmfield(model,field);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
