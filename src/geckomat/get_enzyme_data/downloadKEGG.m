function downloadKEGG(keggID, filePath, keggGeneID)
% downloadKEGG  Download a KEGG organism's gene data and write a CSV.
%
% Downloads gene information for every gene of a KEGG organism code, and
% writes the result to a CSV file. Used by loadDatabases when its
% expected kegg.tsv is missing; can also be called directly to refresh
% one organism's KEGG snapshot without re-running loadDatabases for
% everything else.
%
% Parameters
% ----------
% keggID : char
%     KEGG organism code, e.g. 'sce' for Saccharomyces cerevisiae.
% filePath : char
%     destination CSV path (overwritten if it exists).
% keggGeneID : char
%     which KEGG gene identifier field becomes the output's gene-name
%     column. 'kegg' (or empty, the default) uses the bare KEGG gene id;
%     any other value is looked up as a cross-reference identifier on
%     each KEGG entry (e.g. 'UniProt').
%
% Examples
% --------
%     downloadKEGG('sce', 'kegg.tsv', 'UniProt');
%
% See also
% --------
% loadDatabases, downloadUniProt

if nargin < 3 || isempty(keggGeneID)
    keggGeneID = 'kegg';
end

%% Download gene information
webOptions = weboptions('Timeout',30);
try
    gene_list = webread(['http://rest.kegg.jp/list/' keggID],webOptions);
catch ME
    switch ME.identifier
        case 'MATLAB:webservices:HTTP400StatusCodeError'
            error(['Unable to download data form KEGG with a potentially invalid ID: ' keggID ])
        otherwise
            rethrow(ME)
    end
end
gene_list = regexpi(gene_list, '[^\n]+','match')';
gene_id   = regexpi(gene_list,['(?<=' keggID ':)\S+'],'match');

% Retrieve information for every gene in the list, 10 genes per query
genesPerQuery = 10;
queries = ceil(numel(gene_id)/genesPerQuery);
keggData  = cell(numel(gene_id),1);
PB = progressReport(queries, ['Downloading KEGG data for organism code ' keggID]);
for i = 1:queries
    PB.count;
    % Download batches of genes
    firstIdx = i*genesPerQuery-(genesPerQuery-1);
    lastIdx  = i*genesPerQuery;
    if lastIdx > numel(gene_id) % Last query has probably less genes
        lastIdx = numel(gene_id);
    end
    url      = ['http://rest.kegg.jp/get/' keggID ':' strjoin([gene_id{firstIdx:lastIdx}],['+' keggID ':'])];

    % KEGG may be temporarily unresponsive, retry with exponential backoff
    maxRetries = 5;
    for attempt = 1:maxRetries
        try
            out = webread(url,webOptions);
            break
        catch ME
            % A malformed or non-existing query will not succeed on retry
            noRetry = ismember(ME.identifier,{'MATLAB:webservices:HTTP400StatusCodeError', ...
                                              'MATLAB:webservices:HTTP404StatusCodeError'});
            if noRetry || attempt == maxRetries
                error(['Unable to download KEGG data (attempt %d of %d).\n' ...
                       'URL: %s\n' ...
                       'Check your internet connection and the status of the KEGG ' ...
                       'REST API, then try again.\nOriginal error: %s'], ...
                       attempt, maxRetries, url, ME.message);
            end
            pause(min(2^attempt,30)) % Wait 2, 4, 8 and 16 seconds before retrying
        end
    end
    outSplit = strsplit(out,['///' 10]); %10 is new line character
    if numel(outSplit) < lastIdx-firstIdx+2
        error('KEGG returns less genes per query') %Reduce genesPerQuery
    end
    keggData(firstIdx:lastIdx) = outSplit(1:end-1);
end

%% Parsing of info to keggDB format
% Entries returned by KEGG are separated by '///' followed by an empty line,
% leaving a stray newline at the start of each entry after splitting. This
% would break the startsWith(...,'ENTRY ') checks below, which detect entries
% where the regular expression did not match and returned the entry unchanged
keggData  = strtrim(keggData);

sequence  = regexprep(keggData,'.*AASEQ\s+\d+\s+([A-Z\s])+?\s+NTSEQ.*','$1');
%No AASEQ -> no protein -> not of interest
noProt    = startsWith(sequence,'ENTRY ');
uni       = regexprep(keggData,'.*UniProt: (\S+?)\s.*','$1');
noUni     = startsWith(uni,'ENTRY ');
uni(noProt | noUni)       = [];
keggData(noProt | noUni) = [];
sequence(noProt | noUni)  = [];
sequence  = regexprep(sequence,'\s+','');
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
