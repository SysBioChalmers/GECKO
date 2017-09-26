# -*- coding: utf-8 -*-

from cobra import Model

from geckopy import data


def test_data():
    assert all(x in data.PROTEIN_PROPERTIES.columns for x in ['mw', 'abundance'])
    assert isinstance(data.COBRA_MODELS['full'], Model)
    assert isinstance(data.COBRA_MODELS['batch'], Model)
