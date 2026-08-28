function downloadUniProt(uniprotID, filePath, geneIDfield, type, reviewed)
% downloadUniProt  Download UniProt protein data for an organism and write a TSV.
%
% Downloads accession, gene name, EC number, mass and sequence for every
% protein matching the given organism identifier, and writes the result
% to a TSV file. Used by loadDatabases when its expected uniprot.tsv is
% missing; can also be called directly to refresh one organism's UniProt
% snapshot without re-running loadDatabases for everything else.
%
% Parameters
% ----------
% uniprotID : char
%     UniProt-side identifier whose meaning depends on type, e.g. an NCBI
%     taxonomy ID such as 559292 for Saccharomyces cerevisiae S288C.
% filePath : char
%     destination TSV path (overwritten if it exists).
% geneIDfield : char
%     UniProt field that becomes the second column (the gene matching
%     key) (optional, default 'gene_oln', Ordered Locus Names, which
%     works for most prokaryotes and budding yeast).
% type : char
%     UniProt query field for uniprotID (optional, default
%     'taxonomy_id'). 'taxonomy' is accepted as an alias.
% reviewed : logical
%     whether to restrict the query to reviewed (Swiss-Prot) entries
%     (optional, default true).
%
% Examples
% --------
%     downloadUniProt(559292, 'uniprot.tsv');
%     downloadUniProt(559292, 'uniprot.tsv', 'gene_oln', 'taxonomy_id', true);
%
% See also
% --------
% loadDatabases, downloadKEGG

if nargin < 3 || isempty(geneIDfield)
    geneIDfield = 'gene_oln';
end
if nargin < 4 || isempty(type)
    type = 'taxonomy_id';
end
if nargin < 5 || isempty(reviewed)
    reviewed = true;
end
if strcmp(type,'taxonomy')
    type = 'taxonomy_id';
end
if reviewed
    uniprotRev = 'reviewed:true+AND+';
else
    uniprotRev = '';
end

disp(['Downloading Uniprot data for ' type ' ' num2str(uniprotID) '. This can take a few minutes.'])
url = ['https://rest.uniprot.org/uniprotkb/stream?query=' uniprotRev ...
       type ':' num2str(uniprotID) '&fields=accession%2C' geneIDfield ...
    '%2Cec%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
try
    urlwrite(url,filePath,'Timeout',30);
    fprintf('Model-specific UniProt database stored at %s\n',filePath);
catch
    error(['Download failed, check your internet connection and try again, or manually download: ' url ...
        ' After downloading, store the file as ' filePath])
end
end
