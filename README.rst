.. image:: GECKO.png
   :align: center

|Build Status| |PyPI| |Gitter|

About GECKO
-----------

The **GECKO** toolbox is a Matlab/Python package for enhancing a **G**\ enome-scale model to account for **E**\ nzyme **C**\ onstraints, using **K**\ inetics and **O**\ mics. It is the companion software to the publication:

Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee, Eduard J. Kerkhoven, Jens Nielsen (2017). *Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic constraints.* `Molecular Systems Biology, 13(8):935 <http://www.dx.doi.org/10.15252/msb.20167411>`_

The software comes in two flavors, Python and Matlab scripts to fetch online data and build the published ecYeast7 GECKO models, and a Python package which can be used with `cobrapy <https://opencobra.github.io/cobrapy/>`_ to obtain a ecYeast7 model object, optionally adjusted for provided proteomics data.

Last update: 2018-03-21

This repository is administered by Benjamin J. Sanchez (`@BenjaSanchez <https://github.com/benjasanchez>`_), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.


Building a GECKO model
----------------------

Required software - Python module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `Python 2.7 <https://www.python.org/>`_
- `setuptools for python 2.7 <http://www.lfd.uci.edu/~gohlke/pythonlibs/#setuptools>`_
- SOAPpy:

::

   easy_install-2.7 SOAPpy

Required software - Matlab module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MATLAB <http://www.mathworks.com/>`_ (7.5 or higher) + Optimization Toolbox.
- The `COBRA toolbox for MATLAB <https://github.com/opencobra/cobratoolbox>`_. Note that `libSBML <http://sbml.org/Software/libSBML>`_ and the `SBML toolbox <http://sbml.org/Software/SBMLToolbox>`_ should both be installed. Both of them are free of charge for academic users. Aditionally, you should add the cobra folder to your MATLAB search path.

Usage
~~~~~

See the supporting information of `Sanchez et al. (2017) <https://dx.doi.org/10.15252/msb.20167411>`_


Integrating proteomic data to the yeast GECKO model
---------------------------------------------------

If all you need is the ecYeast7 model to use together with cobrapy you can use the ``geckopy`` Python package.

Required software
~~~~~~~~~~~~~~~~~

- Python 2.7, 3.4, 3.5 or 3.6
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

Contributors
------------

- Benjamin J. Sanchez (`@BenjaSanchez <https://github.com/benjasanchez>`_), Chalmers University of Technology, Gothenburg Sweden
- Ivan Domenzain (`@IVANDOMENZAIN <https://github.com/IVANDOMENZAIN>`_), Chalmers University of Technology, Gothenburg Sweden
- Moritz Emanuel Beber (`@Midnighter <https://github.com/Midnighter>`_), Danish Technical University, Lyngby Denmark
- Henning Redestig (`@hredestig <https://github.com/hredestig>`_), Danish Technical University, Lyngby Denmark
- Cheng Zhang, Science for Life Laboratory, KTH - Royal Institute of Technology, Stockholm Sweden

.. |Build Status| image:: https://travis-ci.org/SysBioChalmers/GECKO.svg?branch=master
   :target: https://travis-ci.org/SysBioChalmers/GECKO
.. |PyPI| image:: https://badge.fury.io/py/geckopy.svg
   :target: https://badge.fury.io/py/geckopy
.. |Gitter| image:: https://badges.gitter.im/SysBioChalmers/GECKO.svg
   :alt: Join the chat at https://gitter.im/SysBioChalmers/GECKO
   :target: https://gitter.im/SysBioChalmers/GECKO?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge