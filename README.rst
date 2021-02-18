.. image:: GECKO.png
   :align: center

|Current Version| |Build Status| |PyPI Version| |Docs Status| |Gitter|

About GECKO
-----------

The **GECKO** toolbox is a Matlab/Python package for enhancing a **G**\ enome-scale model to account for **E**\ nzyme **C**\ onstraints, using **K**\ inetics and **O**\ mics. It is the companion software to `this <http://www.dx.doi.org/10.15252/msb.20167411>`_ publication, and it has two main parts:

- ``geckomat``: Matlab+Python scripts to fetch online data and build/simulate enzyme-constrained models.
- ``geckopy``: a Python package which can be used with `cobrapy <https://opencobra.github.io/cobrapy/>`_ to obtain a ecYeastGEM model object, optionally adjusted for provided proteomics data.

Last update: 2021-02-17

This repository is administered by Benjamin J. Sanchez (`@BenjaSanchez <https://github.com/benjasanchez>`_), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.


geckomat: Building enzyme-constrained models
--------------------------------------------

Required software - Python module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `Python 2.7 <https://www.python.org/>`_
- `setuptools for python 2.7 <http://www.lfd.uci.edu/~gohlke/pythonlibs/#setuptools>`_
- SOAPpy:

  ::

     easy_install-2.7 SOAPpy

Required software - Matlab module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MATLAB <http://www.mathworks.com/>`_ 9.1 (R2016b) or higher + Optimization Toolbox.
- The `COBRA toolbox for MATLAB <https://github.com/opencobra/cobratoolbox>`_.
- The `RAVEN toolbox for MATLAB <https://github.com/SysBioChalmers/RAVEN>`_.
- The `libSBML MATLAB API <https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface>`_ (version 5.17.0 is recommended).

Usage
~~~~~

- **For creating an enzyme constrained model:**

  - Update the following data files in ``/databases`` with your organism infomation:

    - ``databases/prot_abundance.txt``: Protein abundance Data from Pax-DB. If data is not available for your organism, then a relative proteomics dataset (in molar fractions) can be used instead. The required format is a tab-separated file, named as ``databases/relative_proteomics.txt`` , with a single header line and 2 columns; the first with gene IDs and the second with the relative abundances for each protein.
    - ``databases/uniprot.tab``: Gene-proteins data from uniprot.
    - ``databases/chemostatData.tsv``: Chemostat data for estimating GAM (optional, called by ``fitGAM.m``).
    - ``databases/manual_data.txt``: Kcat data from eventual manual curations (optional, called by ``manualModifications.m``).

  - Adapt the following functions in ``/geckomat`` to your organism:

    - ``geckomat/getModelParameters.m``
    - ``geckomat/change_model/manualModifications.m``
    - ``geckomat/limit_proteins/sumProtein.m``
    - ``geckomat/limit_proteins/scaleBioMass.m``
    - ``geckomat/kcat_sensitivity_analysis/changeMedia_batch.m`` (optional)
    - ``geckomat/change_model/removeIncorrectPathways.m`` (optional, called by ``manualModifications.m``)
    - ``geckomat/limit_proteins/sumBioMass.m`` (optional, called by ``sumProtein.m`` & ``scaleBiomass.m``)

  - Run ``geckomat/get_enzyme_data/updateDatabases.m`` to update ``ProtDatabase.mat``.
  - Run ``geckomat/enhanceGEM.m`` with your metabolic model as input.

- **For performing simulations with an enzyme-constrained model:** Enzyme-constrained models can be used as any other metabolic model, with toolboxes such as COBRA or RAVEN. For more information on rxn/met naming convention, see the supporting information of `Sanchez et al. (2017) <https://dx.doi.org/10.15252/msb.20167411>`_

geckopy: Integrating proteomic data to ecYeastGEM
-------------------------------------------------

If all you need is the ecYeastGEM model to use together with cobrapy you can use the ``geckopy`` Python package.

Required software
~~~~~~~~~~~~~~~~~

- Python 3.6, 3.7 or 3.8
- cobrapy

Installation
~~~~~~~~~~~~

::

   pip install geckopy

Usage
~~~~~

.. code:: python

   from geckopy import GeckoModel
   import pandas
   some_measurements = pandas.Series({'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1})
   model = GeckoModel('multi-pool')
   model.limit_proteins(some_measurements)
   model.optimize()

Contributing
------------

Contributions are always welcome! Please read the `contributing guidelines <https://github.com/SysBioChalmers/GECKO/blob/devel/.github/CONTRIBUTING.md>`_ to get started.

Contributors
------------

- Ivan Domenzain (`@IVANDOMENZAIN <https://github.com/IVANDOMENZAIN>`_), Chalmers University of Technology, Gothenburg Sweden
- Eduard Kerkhoven (`@edkerk <https://github.com/edkerk>`_), Chalmers University of Technology, Gothenburg Sweden
- Benjamin J. Sanchez (`@BenjaSanchez <https://github.com/benjasanchez>`_), Chalmers University of Technology, Gothenburg Sweden
- Moritz Emanuel Beber (`@Midnighter <https://github.com/Midnighter>`_), Danish Technical University, Lyngby Denmark
- Henning Redestig (`@hredestig <https://github.com/hredestig>`_), Danish Technical University, Lyngby Denmark
- Cheng Zhang, Science for Life Laboratory, KTH - Royal Institute of Technology, Stockholm Sweden

.. |Current Version| image:: https://badge.fury.io/gh/sysbiochalmers%2Fgecko.svg
   :target: https://badge.fury.io/gh/sysbiochalmers%2Fgecko
.. |Build Status| image:: https://travis-ci.com/SysBioChalmers/GECKO.svg?branch=master
   :target: https://travis-ci.com/SysBioChalmers/GECKO
.. |PyPI Version| image:: https://badge.fury.io/py/geckopy.svg
   :target: https://badge.fury.io/py/geckopy
.. |Docs Status| image:: https://readthedocs.org/projects/geckotoolbox/badge/?version=latest
   :alt: Documentation Status
   :target: http://geckotoolbox.readthedocs.io/
.. |Gitter| image:: https://badges.gitter.im/SysBioChalmers/GECKO.svg
   :alt: Join the chat at https://gitter.im/SysBioChalmers/GECKO
   :target: https://gitter.im/SysBioChalmers/GECKO?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
