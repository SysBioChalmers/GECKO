function downloadKEGG(keggID, filePath, keggGeneID)
% downloadKEGG
%   Download organism-specific protein information from the KEGG REST
%   API and write it to filePath as a comma-delimited table with the
%   schema consumed by loadDatabases.m:
%       uniprot, gene_name, kegg_gene, ec, mw, pathway, sequence
%
%   Extracted from loadDatabases.m (was a local function); now also
%   callable directly to refresh kegg.tsv without going through
%   loadDatabases.
%
% Input:
%   keggID       KEGG organism code (e.g. 'sce')
%   filePath     destination CSV path
%   keggGeneID   KEGG entry field to use as the gene matching key.
%                'kegg' (or empty) copies the bare KEGG gene ID from
%                ENTRY; any other value names an annotation line
%                (e.g. 'OrderedLocus').
%
% Usage:
%   downloadKEGG('sce', '/path/to/data/kegg.tsv', 'kegg')

if nargin < 3 || isempty(keggGeneID)
    keggGeneID = 'kegg';
end

%% Download gene information
progressbar(['Downloading KEGG data for organism code ' keggID])
webOptions = weboptions('Timeout',30);
try
    gene_list = webread(['http://rest.kegg.jp/list/' keggID],webOptions);
catch ME
    switch ME.identifier
        case 'MATLAB:webservices:HTTP400StatusCodeError'
            error(['Unable to download data form KEGG with a potentially invalid ID: ' keggID ])
    end
    rethrow(ME);
end
gene_list = regexpi(gene_list, '[^\n]+','match')';
gene_id   = regexpi(gene_list,['(?<=' keggID ':)\S+'],'match');

% Retrieve information for every gene in the list, 10 genes per query
genesPerQuery = 10;
queries = ceil(numel(gene_id)/genesPerQuery);
keggData  = cell(numel(gene_id),1);
for i = 1:queries
    % Download batches of genes
    firstIdx = i*genesPerQuery-(genesPerQuery-1);
    lastIdx  = i*genesPerQuery;
    if lastIdx > numel(gene_id) % Last query has probably less genes
        lastIdx = numel(gene_id);
    end
    url      = ['http://rest.kegg.jp/get/' keggID ':' strjoin([gene_id{firstIdx:lastIdx}],['+' keggID ':'])];

    retry = true;
    while retry
        try
            retry = false;
            out   = webread(url,webOptions);
        catch
            retry = true;
        end
    end
    outSplit = strsplit(out,['///' 10]); %10 is new line character
    if numel(outSplit) < lastIdx-firstIdx+2
        error('KEGG returns less genes per query') %Reduce genesPerQuery
    end
    keggData(firstIdx:lastIdx) = outSplit(1:end-1);
    progressbar(i/queries)
end

%% Parsing of info to keggDB format
sequence  = regexprep(keggData,'.*AASEQ\s+\d+\s+([A-Z\s])+?\s+NTSEQ.*','$1');
%No AASEQ -> no protein -> not of interest
noProt    = startsWith(sequence,'ENTRY ');
uni       = regexprep(keggData,'.*UniProt: (\S+?)\s.*','$1');
noUni     = startsWith(uni,'ENTRY ');
uni(noProt | noUni)       = [];
keggData(noProt | noUni) = [];
sequence(noProt | noUni)  = [];
sequence  = regexprep(sequence,'\s*','');
keggGene  = regexprep(keggData,'ENTRY\s+(\S+?)\s.+','$1');

switch keggGeneID
    case {'kegg',''}
        gene_name = keggGene;
    otherwise
        % In case there are special characters:
        keggGeneIDT = regexptranslate('escape',keggGeneID);
        gene_name = regexprep(keggData,['.+' keggGeneIDT ': (\S+?)\n.+'],'$1');
        noID = ~contains(keggData,keggGeneIDT);
        if all(noID)
            error(['None of the KEGG entries are annotated with the gene identifier ' keggGeneID])
        else
            gene_name(noID)= [];
            keggData(noID) = [];
            keggGene(noID) = [];
            sequence(noID) = [];
            uni(noID)      = [];
        end
end

EC_names  = regexprep(keggData,'.*ORTHOLOGY.*\[EC:(.*?)\].*','$1');
EC_names(startsWith(EC_names,'ENTRY ')) = {''};

MW = cell(numel(sequence),1);
for i=1:numel(sequence)
    if ~isempty(sequence{i})
        MW{i} = num2str(round(calculateMW(sequence{i})));
    end
end

pathway   = regexprep(keggData,'.*PATHWAY\s+(.*?)(BRITE|MODULE).*','$1');
pathway(startsWith(pathway,'ENTRY ')) = {''};
pathway   = strrep(pathway,[keggID '01100  Metabolic pathways'],'');
pathway   = regexprep(pathway,'\n','');
pathway   = regexprep(pathway,'           ','');

out = [uni, gene_name, keggGene, EC_names, MW, pathway, sequence];
out = cell2table(out);

writetable(out, filePath, 'FileType', 'text', 'WriteVariableNames',false);
fprintf('Model-specific KEGG database stored at %s\n',filePath);
end
