# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os
import pandas as pd
from cobra.io import read_sbml_model


class ModelList(object):
    def __getitem__(self, item):
        key = dict(batch='ecYeast7_batch', full='ecYeast7')[item]
        file_name = os.path.join(os.path.dirname(__file__), 'data_files/{}.xml'.format(key))
        return read_sbml_model(file_name)


"""Should have, for all proteins in model
- uniprot id
- pax abundance in ppm
- molecular weight (or average molecular weight)
"""
ENZYME_PROPERTIES = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data_files/enzymes.txt'), index_col=0)
COBRA_MODELS = ModelList()
