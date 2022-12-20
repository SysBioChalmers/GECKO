from urllib.request import urlretrieve as retrieve
import pandas as pd
import os
from tempfile import mkdtemp

DATADIR = os.path.join(os.path.dirname(__file__), '../geckopy/data_files')


def protein_properties():
    tmpdir = mkdtemp()
    retrieve('http://pax-db.org/data/abundances/4932-WHOLE_ORGANISM-integrated.txt',
             os.path.join(tmpdir, 'yeast-paxdb.txt'))
    retrieve(('http://www.uniprot.org/uniprot/?query='
              'taxonomy%3A559292%20AND%20reviewed%3Ayes&columns=id%2Cmass%2Cgenes(OLN)&format=tab'),
             os.path.join(tmpdir, 'swissprot.txt'))
    pax = pd.read_csv(os.path.join(tmpdir, 'yeast-paxdb.txt'), comment='#', sep='\t',
                      names=['paxid', 'locus', 'abundance'])
    pax = pax[['locus', 'abundance']]
    pax['locus'] = pax['locus'].apply(lambda x: x.replace('4932.', ''))
    pax['abundance'] = pax['abundance'] / 1e6

    mw = pd.read_csv(os.path.join(tmpdir, 'swissprot.txt'), sep='\t', thousands=',')
    mw.columns = ['uniprot', 'mw', 'locus']

    proteins = pd.merge(mw, pax, on='locus', how='outer')
    proteins.index = proteins.uniprot
    (proteins
     [(proteins['uniprot'].notnull())]
     [['mw', 'abundance']]
     .fillna(proteins.mean())
     .to_csv(os.path.join(DATADIR, 'proteins.txt')))

if '__main__' in __name__:
    protein_properties()
