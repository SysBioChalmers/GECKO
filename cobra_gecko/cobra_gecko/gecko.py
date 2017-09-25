# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd
import re
import numpy as np
from six import iteritems
from itertools import chain

from cobra import Reaction, Metabolite, Model

from cobra_gecko.data import ENZYME_PROPERTIES, COBRA_MODELS


class GeckoModel(Model):
    def __init__(self, model, enzyme_properties=None, p_total=0.4005, p_base=0.4005, f=0.4461, sigma=0.5,
                 c_base=0.4067, biomass_reaction='r_4041', protein_pool_exchange='prot_pool_exchange',
                 common_protein_pool='prot_pool'):
        """Convenience class for adjusting gecko-cobra models.

        Parameters
        ----------
        model : cobra.Model
            A cobra model to apply enzyme constraints to.
        enzyme_properties : pd.DataFrame
            A data frame that defined molecular weight (g/mol) 'mw', for 'uniprot' proteins and their average
            'abundance' in ppm.
        p_total : float
            total protein fraction in cell in g protein / g DW. Can be set to a different value than
            p_total to take the higher fraction of ribosomal proteins at higher growth rates into account.
        p_base : float
            protein content at dilution rate 0.1 / h in g protein / g DW.
        f : float
            The fraction of measured proteins versus total proteins in genome (p_model / p_total) (g / g)
        sigma : float
            The parameter adjusting how much of a protein pool can take part in reactions.
        c_base : float
            The carbohydrate content at dilution rate 0.1 / h
        biomass_reaction : str
            The identifier for the biomass reaction
        protein_pool_exchange : str
            The identifier of the protein pool exchange reaction
        common_protein_pool : str
            The identifier of the metabolite representing the common protein pool
        """
        super(GeckoModel, self).__init__(id_or_model=model.copy(), name=model.name)
        self.biomass_reaction = self.reactions.get_by_id(biomass_reaction)
        self.enzyme_properties = enzyme_properties or ENZYME_PROPERTIES
        try:
            self.common_protein_pool = self.metabolites.get_by_id(common_protein_pool)
        except KeyError:
            self.common_protein_pool = Metabolite(common_protein_pool)
        try:
            self.protein_pool_exchange = self.reactions.get_by_id(protein_pool_exchange)
        except KeyError:
            self.protein_pool_exchange = Reaction(protein_pool_exchange)
            self.protein_pool_exchange.add_metabolites({self.common_protein_pool: 1.})
            self.add_reactions([self.protein_pool_exchange])
        self.protein_exchange_re = re.compile(r'^prot_(.*)_exchange$')
        self.pool_protein_exchange_re = re.compile(r'^draw_prot_(.*)$')
        self.concentrations = pd.Series(np.nan, index=self.enzymes)
        self.p_total = p_total
        self.c_base = c_base
        self.p_base = p_base
        self.f_mass_fraction_measured_matched_to_total = f
        self.sigma_saturation_factor = sigma
        self.fp_fraction_protein = self.p_total / self.p_base
        self.fc_carbohydrate_content = (self.c_base + self.p_base - self.p_total) / self.c_base
        self.fs_matched_adjusted = None
        self.p_measured = None
        self.fn_mass_fraction_unmeasured_matched = None
        self.fm_mass_fraction_matched = None

    def fraction_to_ggdw(self, fraction):
        """Convert protein measurements in mass fraction of total to g protein / g DW

        Parameters
        ----------
        fraction : pd.Series
            Data of protein measurements which are absolute quantitative fractions of the total amount of these
            measured proteins. Normalized to sum == 1.

        Returns
        -------
        pd.Series
            g protein / g DW for the measured proteins
        """
        # measurements should be quantitative fractions of the total measured proteins, normalized to unit-length
        fraction = fraction / fraction.sum()
        fraction_measured = self.enzyme_properties['abundance'][list(fraction.index)].sum()
        p_measured = self.p_total * fraction_measured
        return fraction.apply(lambda x: x * p_measured)

    def apply_measurements(self, measurements):
        """Apply proteomics measurements to model.

        Parameters
        ----------
        measurements : pd.Series
            Protein abundances in fraction of total (normalized to sum to 1)
        """
        ggdw = self.fraction_to_ggdw(measurements)
        # * section 2.5
        # 1. define mmmol_gdw as ub for measured enzymes
        for enzyme_id, value in iteritems(ggdw):
            try:
                mmol_gdw = value / (self.enzyme_properties.loc[enzyme_id, 'mw'] / 1000)
                rxn = self.reactions.get_by_id('prot_{}_exchange'.format(enzyme_id))
            except KeyError:
                pass
            else:
                self.concentrations[enzyme_id] = value
                rxn.bounds = 0, mmol_gdw
        # 2. p_measured is aggregate mass of al matched enzymes
        self.p_measured = self.concentrations.sum()
        # 3. fm, mass fraction of measured proteins in the model over total
        self.fm_mass_fraction_matched = self.p_measured / self.p_total
        # 4. mass fraction of unmeasured proteins in the model over all proteins not matched to model
        properties_unmeasured = self.enzyme_properties.loc[self.unmeasured]
        self.fn_mass_fraction_unmeasured_matched = (
            (properties_unmeasured['abundance'] * properties_unmeasured['mw']).sum() /
            (self.enzyme_properties['abundance'] * self.enzyme_properties['mw']).sum())
        self.f_mass_fraction_measured_matched_to_total = self.fn_mass_fraction_unmeasured_matched / (
            1 - self.fm_mass_fraction_matched)
        # 5. constrain unmeasured proteins by common pool
        self.constrain_pool()
        self.adjust_biomass_composition()

    def constrain_pool(self):
        """Constrain common protein pool

        Proteins without their own protein pool are collectively constrained by the common protein pool. Remove
        protein pools for all proteins that don't have measurements, along with corresponding draw reactions,
        and add these to the common protein pool and reaction.
        """
        new_reactions = []
        to_remove = []
        # * section 2.5.1
        # 1. and 2. introduce `prot_pool` and exchange reaction done in __init__
        # 3. limiting total usage with the unmeasured amount of protein
        self.fs_matched_adjusted = ((self.p_total - self.p_measured) * self.f_mass_fraction_measured_matched_to_total *
                                    self.sigma_saturation_factor)
        self.reactions.prot_pool_exchange.bounds = 0, self.fs_matched_adjusted
        # 4. Remove other enzyme usage reactions and replace with pool exchange reactions
        for enzyme_id in self.unmeasured:
            to_remove.extend(self.reactions.query('prot_{}_exchange'.format(enzyme_id)))
            draw_reaction_id = 'draw_prot_{}'.format(enzyme_id)
            if draw_reaction_id not in self.reactions:
                draw_rxn = Reaction(draw_reaction_id)
                protein_pool = self.metabolites.get_by_id('prot_{}_c'.format(enzyme_id))
                metabolites = {self.common_protein_pool: -self.enzyme_properties.loc[enzyme_id, 'mw'] / 1000.,
                               protein_pool: 1}
                draw_rxn.add_metabolites(metabolites)
                new_reactions.append(draw_rxn)
        self.add_reactions(new_reactions)
        self.remove_reactions(to_remove)

    def adjust_biomass_composition(self, gam=31.):
        """Adjust the biomass composition.

        After changing the protein and carbohydrate content based on measurements, adjust the corresponding
        coefficients of the biomass reaction.
        """
        for met in self.biomass_reaction.metabolites:
            coefficient = self.biomass_reaction.metabolites[met]
            sign = -1 if coefficient < 0 else 1
            is_aa = 'tRNA' in met.name
            is_ch = any(x in met.name for x in {'(1->3)-beta-D-glucan', '(1->6)-beta-D-glucan',
                                                'chitin', 'glycogen', 'mannan', 'trehalose'})
            is_atp = 'ATP' in met.name
            is_adp = 'ADP' in met.name
            is_h2o = 'H2O' in met.name
            is_h = 'H+' in met.name
            is_p = 'phosphate' in met.name

            if is_atp or is_adp or is_h2o or is_h or is_p:
                coefficient = sign * (gam + 16.965 * self.fp_fraction_protein + 5.210 * self.fc_carbohydrate_content)
            elif is_aa:
                coefficient = self.fp_fraction_protein * coefficient
            elif is_ch:
                coefficient = self.fc_carbohydrate_content * coefficient
            self.biomass_reaction.metabolites[met] = coefficient

    @property
    def unmeasured(self):
        """Unmeasured enzymes

        Returns
        -------
        list
            The unmeasured enzymes, protein identifiers.
        """
        return list(self.concentrations[self.concentrations.isnull()].index)

    @property
    def enzymes(self):
        return self.individual_enzymes.union(self.pool_enzymes)

    @property
    def individual_enzymes(self):
        """Enzymes with their individual abundance pool.

        Returns
        -------
        frozenset
            The set of proteins that have a defined separate pool exchange reaction.
        """

        return frozenset(chain.from_iterable(re.findall(self.protein_exchange_re, rxn.id)
                                             for rxn in self.protein_exchanges))

    @property
    def pool_enzymes(self):
        """Enzymes

        Returns
        -------
        frozenset
            The set of proteins that have a defined draw reaction.
        """

        return frozenset(chain.from_iterable(re.findall(self.pool_protein_exchange_re, rxn.id)
                                             for rxn in self.protein_exchanges))

    @property
    def protein_exchanges(self):
        """Protein-exchange reactions.

        Returns
        -------
        frozenset
            Set of protein exchange reactions (individual and common protein pool reactions)
        """
        return (frozenset(rxn for rxn in self.reactions
                          if (re.match(self.protein_exchange_re, rxn.id) or
                              re.match(self.pool_protein_exchange_re, rxn.id))) -
                {self.protein_pool_exchange})


def gecko_model(model=None, protein_measurements=None, enzyme_properties=None, p_total=0.4005, p_base=0.4005, f=0.4461,
                sigma=0.5, c_base=0.4067, biomass_reaction='r_4041', protein_pool_exchange='prot_pool_exchange',
                common_protein_pool='prot_pool'):
    """Enzyme constrained metabolic model of yeast.

    Get and adjust an instance of the ecYeast7 model and adjust based on proteomics data.

    Parameters
    ----------
    model : cobra.Model
        A cobra model to apply enzyme constraints to.
    protein_measurements : pd.Series
        A series with protein measurements in fraction of total.
    enzyme_properties : pd.DataFrame
        A data frame that defined molecular weight (g/mol) 'mw', for 'uniprot' proteins and their average
        'abundance' in ppm.
    p_total : float
        total protein fraction in cell in g/gDW
    p_base : float
        protein content at dilution rate 0.1 / h in g/gDW
    f : float
        The fraction of measured proteins versus total proteins in genome (p_model / p_total) (g / g)
    sigma : float
        The parameter adjusting how much of a protein pool can take part in reactions.
    c_base : float
        The carbohydrate content at dilution rate 0.1 / h
    biomass_reaction : str
        The identifier for the biomass reaction
    protein_pool_exchange : str
        The identifier of the protein pool exchange reaction
    common_protein_pool : str
        The identifier of the metabolite representing the common protein pool
    """
    if model is None and protein_measurements is None:
        model = COBRA_MODELS['batch']
    elif model is None:
        model = COBRA_MODELS['full']
    model = GeckoModel(model, enzyme_properties, p_total, p_base, f, sigma, c_base, biomass_reaction,
                       protein_pool_exchange, common_protein_pool)
    if protein_measurements is not None:
        model.apply_measurements(protein_measurements)
    return model
