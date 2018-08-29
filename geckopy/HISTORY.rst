History
=======

1.3.1 (2018-08-28)
------------------
* Features:
    * Adapted the pipeline to work with `yeast-GEM <https://github.com/SysBioChalmers/yeast-GEM>`_, including loading, processing and saving the model. Current model is constructed from yeast `v8.1.3 <https://github.com/SysBioChalmers/yeast-GEM/releases/tag/v8.1.3>`_ (PR #39).
    * When constructing ``ecModel_batch``, lipid fraction is now scaled together with protein and carbohydrate fractions (PR #39).
* Fixes:
    * ``geckopy`` tests flexibilized to comply with yeast-GEM (PR #39).
* Refactoring:
    * Reorganized the repo, making a division between ``geckomat`` (Matlab part for generation + simulation of ecModels) and ``geckopy`` (Python part for simulations of ecYeastGEM) (PR #40).
    * Parameters ``f`` (mass fraction of enzymes in model), ``Pbase``, ``Cbase``, ``Lbase`` (biomass composition) and ``GAM`` (growth-associated ATP maintenance) are now automatically computed (PR #39).
    * Added `RAVEN <https://github.com/SysBioChalmers/RAVEN>`_ as a dependency for ``geckomat`` (PR #38).
    * Changed most COBRA functions in pipeline to RAVEN functions (PR #39).

1.3.0 (2018-08-01)
------------------
* Features:
    * Protein flexibilization: When proteomic measurements are provided, individual protein levels will now be iteratively flexibilized by the pipeline if the model results to be overconstrained, based on a provided growth rate. After this, flexibilized protein exchange pseudoreaction upper bounds will be set to the their flux values from a parsimonious FBA simulation (PR #34).
    * Utilities: Included a folder with useful functions (PR #34).
* Fixes:
    * Fixes #14: CI is no longer failing, as model location, model naming and metabolite ID naming were corrected. ``test_adjust_pool_bounds`` was simplified to test with only 1 essential protein (PR #28).

1.2.1 (2018-05-30)
------------------
* Features:
    * All genes from the original yeast model now included in the ``.xml`` file. Genes connected to enzyme constraints are now stored in ``model.enzGenes`` in the ``.mat`` structure.
    * Docs badge in README.
* Fixes:
    * Fields ``grRules`` and ``rules`` fixed in a consistent way:
        * ``grRules`` for the backwards reactions are the same as for the forward ones.
        * For reactions catalyzed by just 1 enzyme (or complex), ``grRules`` of the original reactions are assigned to them.
        *  For reactions catalyzed by more than 1 enzyme (or more than 1 complex), ``grRules`` of the original reactions are assigned to the arm reactions, and the corresponding sub-rules are assigned to the isozyme-controlled reactions.
        * For enzyme exchange reactions, ``grRules`` are assigned as thecorresponding gene ID.
        * The ``rules`` field is set equal to ``grRules`` for providing consistency with different toolboxes.
    * Inter-OS compatibility:
      * Numbers in scientific notation are stored in the ``.xml`` files with format ``Xe-0N``, not ``Xe-00N``, or with format ``Xe-1N``, not ``Xe-01N``, regardless of the OS used for generating them.
      * Numbers in all files are shown with up to 6 significant figures.
* Refactoring:
    * Updated to new COBRA standards for ``addReaction`` usage.
* NOTE: Not available in pypi (issue #14 unresolved)

1.2.0 (2018-04-12)
------------------
* Implemented automatic *kcat* flexibilization for over-constrained models:
    * Based on a maximum growth rate specified by the user, the algorithm iteratively identifies the top growth-limiting *kcat* value and changes it for the highest one in BRENDA (same EC number)
    * Once that the model is growing close to the set value, the average enzyme saturation factor is refitted
    * For non-feasible/zero-growth models, sensitivity analysis is performed on a reaction and enzyme basis rather than on individual *kcat* values
    * The outputs of this step are stored in ``topUsedEnzymes.txt`` and ``kcatModification.txt`` and can be used for further manual curation
* All databases updated (BRENDA, swissprot, KEGG, PaxDB)
* More generic gene/protein matching for compatibility with other models
* Re-organization of all output files in a single folder
* New badges + styling of website
* NOTE: Not available in pypi (issue #14 unresolved)

1.1.2 (2018-03-20)
------------------
* Improved kcat matching to BRENDA with:
    1) Specific activity
    2) Phylogenetic distance, when data for organism of choice is not available
* Switched to readthedocs for documentation: http://geckotoolbox.readthedocs.io
* Added a Gitter room for discussion: https://gitter.im/SysBioChalmers/GECKO
* Switched to a simplified GitFlow structure (``master`` + ``devel`` + feature branches)
* Python 3.4 environment dropped in CI (no longer supported by pandas)
* NOTE: Not available in pypi (issue #14 unresolved)

1.1.1 (2017-12-08)
------------------
* Model and data are now also deployed.
* Changes in license and readme.

1.1.0 (2017-09-07)
------------------
* First release on PyPI.

1.0.0 (2017-09-07)
------------------
* First release of GECKO in Github.