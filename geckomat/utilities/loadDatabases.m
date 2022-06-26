function [uniprotDB,keggDB] = loadDatabases(taxonID,keggID)

geckoPath = findGECKOroot();
if nargin<2
    param=getModelParameters();
    keggID=param.keggID;
end
if nargin<1
    taxonID=num2str(param.taxonID); % In case it is numeric
end

%Uniprot
if ~exist(fullfile(geckoPath,'databases','uniprotDB.tsv'))
    disp(['Downloading Uniprot data for taxonomic ID ' taxonID '.'])
    url = ['https://www.uniprot.org/uniprot/?query=taxonomy:' taxonID ...
        '&format=tab&compress=no&columns=id,genes(OLN),genes(PREFERRED),ec,database(GeneID),database(RefSeq),mass,sequence'];
    uniprotDL = webread(url);
    fid=fopen(fullfile(geckoPath,'databases','uniprotDB.tsv'),'w');
    fwrite(fid,uniprotDL);
    fclose(fid);
    clear uniprotDL % Load file fresh
else
    disp(['Load existing Uniprot database file for taxonomic ID ' taxonID '.'])
end

fid         = fopen(fullfile(geckoPath,'databases','uniprotDB.tsv'),'r');
fileContent = textscan(fid,'%s %s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);
uniprotDB.ID      = fileContent{1};
uniprotDB.genes   = fileContent{2};
uniprotDB.eccodes = fileContent{4};
uniprotDB.MW      = str2double(fileContent{7});
uniprotDB.seq     = fileContent{8};

% KEGG
% keggID     = parameters.keggID;
% cd (current)
% if ~regexp(keggID,'[a-z]{3,4}')
%     error('Please specify the KEGG organism ID in the script getModelParameters.m')
% end
% 
% %Build Swissprot table:
% swissprot = buildSWISSPROTtable;
% 
% %Download KEGG data:
% mkdir ../../databases/KEGG
% downloadKEGGdata(keggID)
% 
% %Build KEGG table
% kegg = buildKEGGtable(keggID);
end

