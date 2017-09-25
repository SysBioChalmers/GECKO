# -*- coding: utf-8 -*-

from cobra_gecko import gecko_model
from cobra_gecko.data import ENZYME_PROPERTIES as props

import pandas as pd

in_model = {'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1, 'P39708': 0.1, 'P39714': 0.1, 'P39726': 0.1, 'Q01574': 0.1}
not_in_model = {'P10591': 0.1, 'P31383': 0.1, 'P32471': 0.1}


def test_gecko_adjustment():
    measurements = pd.concat([pd.Series(in_model), pd.Series(not_in_model)])
    model = gecko_model(protein_measurements=measurements)
    sol = model.optimize()
    assert sol.objective_value > 0.01
    # assert all(rxn.ub == props[] for rxn in )
