function [uniprotDB,keggDB] = loadDatabases(taxonID,keggID)

if nargin<1
    %Use organism from GECKOModelAdapter
    global GECKOModelAdapter
    param=checkGECKOModelAdapter;
    keggID=param.keggID;
    taxonID=param.keggID;
end

geckoPath = findGECKOroot();
param=modelAdapter.getParameters();
if nargin<3
    keggID=param.keggID;
end
if nargin<2
    taxonID=num2str(param.taxonID); % In case it is numeric
end

%Uniprot
filename = [num2str(taxonID) '.tsv'];
if ~exist(fullfile(geckoPath,'databases/uniprot/',filename))
    disp(['Downloading Uniprot data for taxonomic ID ' taxonID '. This can take a few minutes.'])
    %url = ['https://www.uniprot.org/uniprot/?query=taxonomy:' taxonID ...
    %    '&format=tab&compress=no&columns=id,genes(OLN),genes(PREFERRED),ec,database(GeneID),database(RefSeq),mass,sequence'];
    url = ['https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:' num2str(taxonID) ...
           '&fields=accession%2C' param.uniprotGeneIdField ...
           '%2Cgene_primary%2Cec%2Cxref_geneid%2Cxref_refseq%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
    urlwrite(url,fullfile(geckoPath,'databases/uniprot/',filename));
%    options=weboptions('Timeout',100); %If the code fails, this may need to be increased
%    uniprotDL = webread(url,options);
%    dl2 = decode(uniprotDL);
%    fid=fopen(fullfile(geckoPath,'databases/uniprot/',filename),'w');
%    fwrite(fid,uniprotDL);
%    fclose(fid);
%    gunzip(fullfile(geckoPath,'databases/uniprot/',filename))
%    clear uniprotDL % Load file fresh
else
    disp(['Loading existing Uniprot database file for taxonomic ID ' taxonID '.'])
end

fid         = fopen(fullfile(geckoPath,'databases/uniprot/',filename),'r');
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
%{
url = 'https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:9606&fields=accession%2Cgene_oln%2Cgene_primary%2Cec%2Cxref_geneid%2Cxref_refseq%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc';
url2 = 'https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:9606&fields=accession%2Cgene_oln%2Cgene_primary%2Cec%2Cxref_geneid%2Cxref_refseq%2Cmass%2Csequence&format=tsv&compressed=true&sort=protein_name%20asc';

tic
urlwrite(url,'C:\Work\MatlabCode\components\Gecko\Gecko3\GECKO\Databases\uniprottest.tsv');
toc % 96 s

tic
urlwrite(url2,'C:\Work\MatlabCode\components\Gecko\Gecko3\GECKO\Databases\uniprottest2.tsv');
toc %115 s
%}