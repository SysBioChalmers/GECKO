# -*- coding: utf-8 -*-

"""Provides data that is shipped with the geckopy package. Copied from the main repository on package build."""

from __future__ import absolute_import

import os
import pandas as pd
from math import isinf
from cobra.io import read_sbml_model


class ModelList(object):
    """List of shipped GECKO models.

    Implements lazy loading, models are only loaded from disk when first requested.
    """

    models = {}

    def __getitem__(self, item):
        """Get a bundled GECKO model.

        Parameters
        ----------
        item : basestring
            Either 'batch' for the single-protein pool ecYeast7 model or 'full' for individually modeled protein pools.

        """
        key = dict(batch='ecYeast7_batch', full='ecYeast7')[item]
        file_name = os.path.join(os.path.dirname(__file__), 'data_files/{}.xml'.format(key))
        model = read_sbml_model(file_name)
        for rxn in model.reactions:
            if rxn.upper_bound == inf:
                rxn.upper_bound = 1000
        return model


"""Should have, for all proteins in model
- uniprot id
- pax abundance in ppm
- molecular weight (or average molecular weight)
"""
PROTEIN_PROPERTIES = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data_files/proteins.txt'), index_col=0)
COBRA_MODELS = ModelList()
