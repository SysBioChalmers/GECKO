function downloadUniProt(uniprotID, filePath, idType, geneIDfield, reviewed)
% downloadUniProt
%   Download organism-specific UniProt data via the REST API and write
%   it to filePath as a tab-delimited TSV with the schema consumed by
%   loadDatabases.m:
%       Entry, <geneIDfield>, EC number, Mass, Sequence
%
%   Extracted from loadDatabases.m so the download step is callable
%   directly (parallels downloadKEGG).
%
% Input:
%   uniprotID    UniProt query identifier (e.g. NCBI taxonomy id 559292
%                when idType is 'taxonomy_id')
%   filePath     destination TSV path
%   idType       UniProt query field for uniprotID. Defaults to
%                'taxonomy_id'. The alias 'taxonomy' is accepted.
%   geneIDfield  UniProt field used as the Gene Names column. Defaults
%                to 'gene_oln'.
%   reviewed     true (default) restricts to reviewed (Swiss-Prot)
%                entries.
%
% Usage:
%   downloadUniProt(559292, '/path/to/data/uniprot.tsv', 'taxonomy_id', 'gene_oln', true)

if nargin < 3 || isempty(idType)
    idType = 'taxonomy_id';
end
if nargin < 4 || isempty(geneIDfield)
    geneIDfield = 'gene_oln';
end
if nargin < 5 || isempty(reviewed)
    reviewed = true;
end
if strcmp(idType,'taxonomy')
    idType = 'taxonomy_id';
end
if reviewed
    uniprotRev = 'reviewed:true+AND+';
else
    uniprotRev = '';
end

disp(['Downloading Uniprot data for ' idType ' ' num2str(uniprotID) '. This can take a few minutes.'])
url = ['https://rest.uniprot.org/uniprotkb/stream?query=' uniprotRev ...
       idType ':' num2str(uniprotID) '&fields=accession%2C' geneIDfield ...
    '%2Cec%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
try
    urlwrite(url,filePath,'Timeout',30);
    fprintf('Model-specific UniProt database stored at %s\n',filePath);
catch
    error(['Download failed, check your internet connection and try again, or manually download: ' url ...
        ' After downloading, store the file as ' filePath])
end
end
