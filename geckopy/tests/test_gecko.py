# -*- coding: utf-8 -*-

from geckopy import GeckoModel

import pandas as pd

in_model = {'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1, 'P39708': 0.1, 'P39714': 0.1, 'P39726': 0.1, 'Q01574': 0.1}
not_in_model = {'P10591': 0.1, 'P31383': 0.1, 'P32471': 0.1}


def test_gecko_adjustment():
    measurements = pd.concat([pd.Series(in_model), pd.Series(not_in_model)])
    model = GeckoModel('multi-pool', protein_measurements=pd.Series(measurements))
    sol = model.optimize()
    assert sol.objective_value > 0.05
    assert len(model.enzymes) - len(model.pool_enzymes) - len(in_model) == 0
    # fraction_measured = model.protein_properties['abundance'][list(fraction.index)].sum()
    # assert all(rxn.ub == props[] for rxn in )
