# -*- coding: utf-8 -*-
"""Implement the GeckoModel which subclasses cobrapy's Model."""
from __future__ import absolute_import

import pandas as pd
import re
import numpy as np
from six import iteritems, string_types
from sympy import S
from itertools import chain

from cobra import Reaction, Metabolite, Model

from geckopy.data import PROTEIN_PROPERTIES, COBRA_MODELS


class GeckoModel(Model):
    """Class for representing GECKO models.

    Implement a model class for Genome-scale model to account for Enzyme Constraints, using Kinetics and Omics [1]_.

    Parameters
    ----------
    model : cobra.Model, str
        A cobra model to apply protein constraints to. Can be 'single-pool' for the bundled ecYeast7 model using only
        a single pool for all proteins, or 'multi-pool' for the model that has separate pools for all measured proteins.
    protein_properties : pd.DataFrame
        A data frame that defined molecular weight (g/mol) 'mw', for 'uniprot' proteins and their average
        'abundance' in ppm.
    sigma : float
        The parameter adjusting how much of a protein pool can take part in reactions. Fitted parameter, default is
        optimized for chemostat experiment in [1]_.
    gam : float
        The growth associated maintenance cost in mmol / gDW. Default fitted for yeast 8.1.3.
    amino_acid_polymerization_cost : float
        The cost for turning amino-acids in proteins in mmol / g. Default taken from [2]_.
    carbohydrate_polymerization_cost : float
        The cost for turning monosaccharides in polysaccharides in mmol / g. Default taken from [2]_.
    c_base : float
        The carbohydrate content at dilution rate 0.1 / h. Default taken from yeast 8.1.3.
    biomass_reaction_id : str
        The identifier for the biomass reaction
    protein_reaction_id : str
        The identifier for the protein reaction
    carbohydrate_reaction_id : str
        The identifier for the carbohydrate reaction
    protein_pool_exchange_id : str
        The identifier of the protein pool exchange reaction
    common_protein_pool_id : str
        The identifier of the metabolite representing the common protein pool

    References
    ----------
    .. [1] Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee, Eduard J. Kerkhoven, Jens Nielsen (
       2017). Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic
       constraints. [Molecular Systems Biology, 13(8): 935, http://www.dx.doi.org/10.15252/msb.20167411

       [2] J. Förster, I. Famili, B. Ø. Palsson and J. Nielsen, Genome Res., 2003, 244–253.

    """

    def __init__(self, model, protein_properties=None,
                 sigma=0.46, c_base=0.3855, gam=36.6, amino_acid_polymerization_cost=37.7,
                 carbohydrate_polymerization_cost=12.8, biomass_reaction_id='r_4041',
                 protein_reaction_id='r_4047', carbohydrate_reaction_id='r_4048',
                 protein_pool_exchange_id='prot_pool_exchange', common_protein_pool_id='prot_pool'):
        """Get a new GECKO model object."""
        model = COBRA_MODELS[model].copy() if isinstance(model, string_types) else model
        super(GeckoModel, self).__init__(id_or_model=model, name=model.name)
        self.biomass_reaction = self.reactions.get_by_id(biomass_reaction_id)
        self.protein_reaction = self.reactions.get_by_id(protein_reaction_id)
        self.carbohydrate_reaction = self.reactions.get_by_id(carbohydrate_reaction_id)
        self.protein_properties = protein_properties or PROTEIN_PROPERTIES
        try:
            self.common_protein_pool = self.metabolites.get_by_id(common_protein_pool_id)
        except KeyError:
            self.common_protein_pool = Metabolite(common_protein_pool_id)
        try:
            self.protein_pool_exchange = self.reactions.get_by_id(protein_pool_exchange_id)
        except KeyError:
            self.protein_pool_exchange = Reaction(protein_pool_exchange_id)
            self.protein_pool_exchange.add_metabolites({self.common_protein_pool: 1.})
            self.add_reactions([self.protein_pool_exchange])
        self.protein_exchange_re = re.compile(r'^prot_(.*)_exchange$')
        self.pool_protein_exchange_re = re.compile(r'^draw_prot_(.*)$')
        self.concentrations = pd.Series(np.nan, index=self.proteins)
        self.gam = gam
        self.amino_acid_polymerization_cost = amino_acid_polymerization_cost
        self.carbohydrate_polymerization_cost = carbohydrate_polymerization_cost
        self.c_base = c_base
        self.sigma_saturation_factor = sigma
        self.measured_ggdw = None
        self.fp_fraction_protein = None
        self.fc_carbohydrate_content = None
        self.p_total = None
        self.c_total = None
        self.p_base = None
        self.fn_mass_fraction_unmeasured_matched = None
        self.fs_matched_adjusted = None
        self.p_measured = None
        self.f_mass_fraction_measured_matched_to_total = None
        self.fm_mass_fraction_matched = None

    def fraction_to_ggdw(self, fraction):
        """Convert protein measurements in mass fraction of total to g protein / g DW.

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
        fraction_measured = self.protein_properties['abundance'][list(fraction.index)].sum()
        p_measured = self.p_total * fraction_measured
        return fraction.apply(lambda x: x * p_measured)

    def limit_proteins(self, fractions=None, ggdw=None, p_total=0.448, p_base=0.46):
        """Apply proteomics measurements to model.

        Apply measurements in the form of fractions of total of the measured proteins, or directly as g / gDW. Must
        supply exactly one of `fractions` or `ggdw`.

        Parameters
        ----------
        fractions : pd.Series
            Protein abundances in fraction of total (normalized to sum to 1). Ignored if `ggdw` is also supplied.
        ggdw : pd.Series
            Protein abundances in g / gDW
        p_total : float
            measured total protein fraction in cell in g protein / g DW. Should be measured for each experiment,
            the default here is taken from [1]_.
        p_base : float
            protein content at dilution rate 0.1 / h in g protein / g DW. Default taken from yeast 8.1.3.

        References
        ----------
        .. [1] Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee, Eduard J. Kerkhoven, Jens Nielsen (
           2017). Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic
           constraints. [Molecular Systems Biology, 13(8): 935, http://www.dx.doi.org/10.15252/msb.20167411


        """
        self.p_total = p_total
        self.p_base = p_base
        self.c_total = self.c_base + self.p_base - self.p_total
        self.fp_fraction_protein = self.p_total / self.p_base
        self.fc_carbohydrate_content = self.c_total / self.c_base
        self.measured_ggdw = self.fraction_to_ggdw(fractions) if ggdw is None else ggdw
        # * section 2.5
        # 1. define mmmol_gdw as ub for measured proteins
        for protein_id, value in iteritems(self.measured_ggdw):
            try:
                mmol_gdw = value / (self.protein_properties.loc[protein_id, 'mw'] / 1000)
                rxn = self.reactions.get_by_id('prot_{}_exchange'.format(protein_id))
                rxn.annotation['uniprot'] = protein_id
            except KeyError:
                pass
            else:
                self.concentrations[protein_id] = value
                rxn.bounds = 0, mmol_gdw
        # 2. p_measured is aggregate mass of all matched proteins
        self.p_measured = self.concentrations.sum()
        # 3. fm, mass fraction of measured proteins in the model over total
        self.fm_mass_fraction_matched = self.p_measured / self.p_total
        # 4. mass fraction of unmeasured proteins in the model over all proteins not matched to model
        self.fn_mass_fraction_unmeasured_matched = (
            self.protein_properties.loc[self.unmeasured_proteins].prod(axis=1).sum() /
            self.protein_properties.prod(axis=1).sum()
        )
        self.f_mass_fraction_measured_matched_to_total = (
            self.fn_mass_fraction_unmeasured_matched / (1 - self.fm_mass_fraction_matched))
        # 5. constrain unmeasured proteins by common pool
        self.constrain_pool()
        self.adjust_biomass_composition()

    def constrain_pool(self):
        """Constrain the draw reactions for the unmeasured (common protein pool) proteins.

        Proteins without their own protein pool are collectively constrained by the common protein pool. Remove
        protein pools for all proteins that don't have measurements, along with corresponding draw reactions,
        and add these to the common protein pool and reaction.
        """
        new_reactions = []
        to_remove = []
        # * section 2.5.1
        # 1. and 2. introduce `prot_pool` and exchange reaction done in __init__
        # 3. limiting total usage with the unmeasured amount of protein
        # looks like the matlab code:
        # self.fs_matched_adjusted = ((self.p_total - self.p_measured) / self.p_base *
        #                             self.f_mass_fraction_measured_matched_to_total *
        #                             self.sigma_saturation_factor)
        # but this gives results more like reported:
        self.fs_matched_adjusted = ((self.p_total - self.p_measured) *
                                    self.f_mass_fraction_measured_matched_to_total *
                                    self.sigma_saturation_factor)
        self.reactions.prot_pool_exchange.bounds = 0, self.fs_matched_adjusted
        # 4. Remove other enzyme usage reactions and replace with pool exchange reactions
        average_mmw = self.protein_properties['mw'].mean() / 1000.
        for protein_id in self.unmeasured_proteins:
            to_remove.extend(self.reactions.query('prot_{}_exchange'.format(protein_id)))
            draw_reaction_id = 'draw_prot_{}'.format(protein_id)
            if draw_reaction_id not in self.reactions:
                draw_rxn = Reaction(draw_reaction_id)
                draw_rxn.annotation['uniprot'] = protein_id
                protein_pool = self.metabolites.get_by_id('prot_{}_c'.format(protein_id))
                try:
                    mmw = self.protein_properties.loc[protein_id, 'mw'] / 1000.
                except KeyError:
                    mmw = average_mmw
                metabolites = {self.common_protein_pool: -mmw, protein_pool: 1}
                draw_rxn.add_metabolites(metabolites)
                new_reactions.append(draw_rxn)
        self.add_reactions(new_reactions)
        self.remove_reactions(to_remove)

    def adjust_biomass_composition(self):
        """Adjust the biomass composition.

        After changing the protein and carbohydrate content based on measurements, adjust the corresponding
        coefficients of the biomass reaction.
        """
        for met in self.protein_reaction.metabolites:
            is_prot = 'protein' in met.name
            if not is_prot:
                coefficient = self.fp_fraction_protein * self.protein_reaction.metabolites[met]
                self.protein_reaction.metabolites[met] = coefficient

        for met in self.carbohydrate_reaction.metabolites:
            is_carb = 'carbohydrate' in met.name
            if not is_carb:
                coefficient = self.fc_carbohydrate_content * self.carbohydrate_reaction.metabolites[met]
                self.carbohydrate_reaction.metabolites[met] = coefficient

        for met in self.biomass_reaction.metabolites:
            sign = -1 if self.biomass_reaction.metabolites[met] < 0 else 1
            is_atp = 'ATP' in met.name
            is_adp = 'ADP' in met.name
            is_h2o = 'H2O' in met.name
            is_h = 'H+' in met.name
            is_p = 'phosphate' in met.name
            if is_atp or is_adp or is_h2o or is_h or is_p:
                coefficient = sign * (self.gam +
                                      self.amino_acid_polymerization_cost * self.p_total +
                                      self.carbohydrate_polymerization_cost * self.c_total)
                self.biomass_reaction.metabolites[met] = coefficient

    def adjust_pool_bounds(self, min_objective=0.05, inplace=False, tolerance=1e-9):
        """Adjust protein pool bounds minimally to make model feasible.

        Bounds from measurements can make the model non-viable or even infeasible. Adjust these minimally by minimizing
        the positive deviation from the measured values.

        Parameters
        ----------
        min_objective : float
            The minimum value of for the ojective for calling the model viable.
        inplace : bool
            Apply the adjustments to the model.
        tolerance : float
            Minimum non-zero value. Solver specific value.

        Returns
        -------
        pd.DataFrame
            Data frame with the series 'original' bounds and the new 'adjusted' bound, and the optimized 'addition'.

        """
        with self as model:
            problem = model.problem
            constraint_objective = problem.Constraint(model.objective.expression, name='constraint_objective',
                                                      lb=min_objective)
            to_add = [constraint_objective]
            new_objective = S.Zero
            for pool in model.individual_protein_exchanges:
                ub_diff = problem.Variable('pool_diff_' + pool.id, lb=0, ub=None)
                current_ub = problem.Variable('measured_bound_' + pool.id, lb=pool.upper_bound, ub=pool.upper_bound)
                constraint = problem.Constraint(pool.forward_variable - current_ub - ub_diff, ub=0,
                                                name='pool_ub_' + pool.id)
                to_add.extend([ub_diff, current_ub, constraint])
                new_objective += ub_diff
                pool.bounds = 0, 1000.
            model.add_cons_vars(to_add)
            model.objective = problem.Objective(new_objective, direction='min')
            model.slim_optimize(error_value=None)
            primal_values = model.solver.primal_values
        adjustments = [(pool.id, primal_values['pool_diff_' + pool.id], pool.upper_bound)
                       for pool in model.individual_protein_exchanges
                       if primal_values['pool_diff_' + pool.id] > tolerance]
        result = pd.DataFrame(adjustments, columns=['reaction', 'addition', 'original'])
        result['adjusted'] = result['addition'] + result['original']
        if inplace:
            for adj in result.itertuples():
                model.reactions.get_by_id(adj.reaction).upper_bound = adj.adjusted
        return result

    @property
    def measured_proteins(self):
        """Get the identifiers of the measured proteins.

        Returns
        -------
        frozenset
            The identifiers for the unmeasured proteins.

        """
        return frozenset(self.concentrations[self.concentrations.notnull()].index)

    @property
    def unmeasured_proteins(self):
        """Get the identifiers of the proteins .

        Returns
        -------
        frozenset
            The protein identifiers for the measured proteins.

        """
        return frozenset(self.concentrations[self.concentrations.isnull()].index)

    @property
    def proteins(self):
        """Get all proteins.

        Returns
        -------
        frozenset
           The set of all proteins identifiers.

        """
        return self.individual_proteins.union(self.pool_proteins)

    @property
    def individual_proteins(self):
        """Get the identifiers for the proteins with their individual abundance pool.

        Returns
        -------
        frozenset
            The set of proteins that have a defined separate pool exchange reaction.

        """
        return frozenset(chain.from_iterable(re.findall(self.protein_exchange_re, rxn.id)
                                             for rxn in self.protein_exchanges))

    @property
    def pool_proteins(self):
        """Get proteins modeled by common protein pool.

        Returns
        -------
        frozenset
            The set of proteins that have a defined draw reaction.

        """
        return frozenset(chain.from_iterable(re.findall(self.pool_protein_exchange_re, rxn.id)
                                             for rxn in self.protein_exchanges))

    @property
    def individual_protein_exchanges(self):
        """Individual protein-exchange reactions.

        Returns
        -------
        frozenset
            Set of protein exchange reactions with individual pools

        """
        return (frozenset(rxn for rxn in self.reactions
                          if re.match(self.protein_exchange_re, rxn.id)) -
                {self.protein_pool_exchange})

    @property
    def pool_protein_exchanges(self):
        """Protein-exchange reactions by single pool.

        Returns
        -------
        frozenset
            Set of protein exchange reactions for single pool reactions.

        """
        return (frozenset(rxn for rxn in self.reactions
                          if re.match(self.pool_protein_exchange_re, rxn.id)) -
                {self.protein_pool_exchange})

    @property
    def protein_exchanges(self):
        """Protein-exchange reactions.

        Returns
        -------
        frozenset
            Set of protein exchange reactions (individual and common protein pool reactions)

        """
        return self.individual_protein_exchanges.union(self.pool_protein_exchanges)
