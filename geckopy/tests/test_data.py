# -*- coding: utf-8 -*-

import pytest
from cobra import Model

from geckopy import data


def test_data():
    assert all(x in data.PROTEIN_PROPERTIES.columns for x in ['mw', 'abundance'])
    assert isinstance(data.COBRA_MODELS['single-pool'], Model)
    assert isinstance(data.COBRA_MODELS['multi-pool'], Model)
    with pytest.raises(KeyError):
        print(data.COBRA_MODELS['blubb'])
