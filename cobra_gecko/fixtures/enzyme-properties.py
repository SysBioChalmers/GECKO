import ssl
from urllib.request import urlopen, urlretrieve, URLError
import pandas as pd
from gzip import GzipFile
import os
from tempfile import mkdtemp


DATADIR = os.path.join(os.path.dirname(__file__), '../cobra_gecko/data_files')


def retrieve(url, f_name):
    try:
        urlretrieve(url, f_name)
    except URLError:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = eval('ssl.CERT_NONE')
        with urlopen(url, context=ctx) as u, open(f_name, 'w') as f:
            f.write(u.read().decode())


def enzyme_properties():
    tmpdir = mkdtemp()
    retrieve(('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping'
              '/by_organism/YEAST_559292_idmapping.dat.gz'),
             os.path.join(tmpdir, 'idmapping.dat.gz'))

    retrieve('http://pax-db.org/data/abundances/4932-WHOLE_ORGANISM-integrated.txt',
             os.path.join(tmpdir, 'yeast-paxdb.txt'))
    retrieve('http://downloads.yeastgenome.org/curation/calculated_protein_info/protein_properties.tab',
             os.path.join(tmpdir, 'protein_properties.tab'))
    idmap = pd.read_csv(GzipFile(os.path.join(tmpdir, 'idmapping.dat.gz'), 'r'), sep='\t',
                        names=['uniprot', 'key', 'locus'])
    idmap = idmap[(idmap.key == 'Gene_OrderedLocusName')][['uniprot', 'locus']]

    # pax = pd.read_csv(os.path.join(tmpdir, 'yeast-paxdb.txt'), comment='#', sep='\t',
    #                   names=['paxid', 'locus', 'abundance'])
    pax = pd.read_csv('/Users/henred/proj/gecko-model/Matlab_Module/limit_proteins/prot_abundance.txt',
                      comment='#', sep='\t',
                      names=['paxid', 'locus', 'abundance'])
    pax = pax[['locus', 'abundance']]
    pax['locus'] = pax['locus'].apply(lambda x: x.replace('4932.', ''))
    pax['abundance'] = pax['abundance'] / 1e6

    mw = pd.read_csv(os.path.join(tmpdir, 'protein_properties.tab'), sep='\t')[['ORF', 'Mw']]
    mw.columns = ['locus', 'mw']

    enzyme = pd.merge(mw, idmap, on='locus')
    enzyme = pd.merge(enzyme, pax, on='locus')
    enzyme = enzyme[['uniprot', 'mw', 'abundance']]
    enzyme.index = enzyme.uniprot
    enzyme = enzyme[['mw', 'abundance']][(enzyme.notnull())]
    enzyme.to_csv(os.path.join(DATADIR, 'enzymes.txt'))


if '__main__' in __name__:
    enzyme_properties()
