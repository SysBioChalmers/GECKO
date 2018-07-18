# -*- coding: utf-8 -*-

"""Provides data that is shipped with the geckopy package. Copied from the main repository on package build."""

from __future__ import absolute_import

import os
import re
import pandas as pd
from math import isinf
from cobra.io import read_sbml_model

DATA_FILES = os.path.join(os.path.dirname(__file__), 'data_files')


class ModelList(object):
    """List of shipped GECKO models.

    Implements lazy loading, models are only loaded from disk when first requested.
    """

    models = {}
    model_files = dict((re.findall(r'_(.*).xml', f)[0], f) for f in os.listdir(DATA_FILES) if f.endswith('.xml'))

    def __getitem__(self, item):
        """Get a bundled GECKO model.

        Parameters
        ----------
        item : basestring
            Either 'single-pool' for the single-protein pool ecYeastGEM model or 'multi-pool' for individually modeled
            protein pools.

        """
        try:
            file_name = self.model_files[item]
        except KeyError:
            raise KeyError('model name must be one of {}'.format(', '.join(list(self.model_files))))
        if file_name not in self.models:
            model = read_sbml_model(os.path.join(os.path.dirname(__file__), 'data_files/{}'.format(file_name)))
            for met in model.metabolites:
                met.id = met.id.replace('__91__', '_')
                met.id = met.id.replace('__93__', '')
            for rxn in model.reactions:
                if isinf(rxn.upper_bound):
                    rxn.upper_bound = 1000
            self.models[file_name] = model
        return self.models[file_name]


"""Should have, for all proteins in model
- uniprot id
- pax abundance in ppm
- molecular weight (or average molecular weight)
"""
PROTEIN_PROPERTIES = pd.read_csv(os.path.join(DATA_FILES, 'proteins.txt'), index_col=0)
COBRA_MODELS = ModelList()
