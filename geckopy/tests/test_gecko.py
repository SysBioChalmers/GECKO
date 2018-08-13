# -*- coding: utf-8 -*-

from geckopy import GeckoModel
from geckopy.data import PROTEIN_PROPERTIES

import os
import pandas as pd


def test_basic_gecko_adjustment():
    in_model = {'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1, 'P39708': 0.1, 'P39714': 0.1, 'P39726': 0.1, 'Q01574': 0.1}
    not_in_model = {'P10591': 0.1, 'P31383': 0.1, 'P32471': 0.1}
    measurements = pd.concat([pd.Series(in_model), pd.Series(not_in_model)])
    model = GeckoModel('multi-pool')
    model.limit_proteins(fractions=pd.Series(measurements))
    sol = model.optimize()
    assert sol.objective_value > 0.05
    assert len(model.proteins) - len(model.pool_proteins) - len(in_model) == 0
    assert all(rxn.upper_bound > 0 for rxn in model.individual_protein_exchanges)


def test_gecko_adjustment_sanchez_etal():
    mmol_gdw = pd.Series.from_csv(os.path.join(os.path.dirname(__file__), '../geckopy/data_files/sanchez-mmol_gdw.csv'))
    ggdw = pd.Series(PROTEIN_PROPERTIES.loc[mmol_gdw.index, 'mw'] / 1000.) * pd.Series(mmol_gdw)
    model = GeckoModel('multi-pool')
    growth_rate_unlimited_protein = model.slim_optimize()
    model.limit_proteins(ggdw=pd.Series(ggdw))
    growth_rate_limited_protein = model.slim_optimize()
    # should be smaller, but how much..
    assert growth_rate_limited_protein < 0.8 * growth_rate_unlimited_protein
    measured_in_model = set(mmol_gdw.index).intersection(model.proteins)
    assert sum(model.concentrations[p] - ggdw[p] for p in measured_in_model) < 1e-10
    assert sum(abs(rxn.upper_bound - mmol_gdw[rxn.annotation['uniprot']])
               for rxn in model.individual_protein_exchanges) < 1e-6
    assert sum(rxn.metabolites[model.common_protein_pool] +
               PROTEIN_PROPERTIES.loc[rxn.annotation['uniprot'], 'mw'] / 1000.
               for rxn in model.pool_protein_exchanges) < 1e-6
    assert model.p_measured > 0.25                                  # With yeast 8.1.3 -> p_measured = 0.296
    assert model.f_mass_fraction_measured_matched_to_total > 0.25   # With yeast 8.1.3 -> f = 0.304
    assert model.protein_pool_exchange.upper_bound > 0.015          # With yeast 8.1.3 -> pool_exchange = 0.0212


def test_adjust_pool_bounds():
    essential = {'P00498': 0.}
    in_model = {'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1, 'P39708': 0.1, 'P39714': 0.1, 'P39726': 0.1, 'Q01574': 0.1}
    expected = set('prot_{}_exchange'.format(pool_id) for pool_id in essential)
    measurements = pd.concat([pd.Series(in_model), pd.Series(essential)])
    model = GeckoModel('multi-pool')
    model.limit_proteins(fractions=pd.Series(measurements))
    assert model.slim_optimize() < 1e-3
    model.adjust_pool_bounds(inplace=False)
    assert model.slim_optimize() < 1e-3
    adjustment = model.adjust_pool_bounds(inplace=True)
    observed = set(adjustment['reaction'])
    assert abs(model.slim_optimize() - 0.05) < 1e-3
    assert expected == observed
