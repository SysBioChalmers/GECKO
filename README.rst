The **GECKO** toolbox is a Matlab/Python package for enhancing a
**G**\ enome-scale model to account for **E**\ nzyme **C**\ onstraints,
using **K**\ inetics and **O**\ mics. It is the companion software to
the following publication:

*Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee,
Eduard J. Kerkhoven, Jens Nielsen (2017). Improving the phenotype
predictions of a yeast genome-scale metabolic model by incorporating
enzymatic constraints. `Molecular Systems Biology, 13(8):
935 <http://www.dx.doi.org/10.15252/msb.20167411>`__ *

GECKO was written by Benjamin J. Sanchez (@BenjaSanchez), Division of
Systems and Synthetic Biology, Department of Biology and Biological
Engineering, Chalmers University of Technology, and Cheng Zhang, Science
for Life Laboratory, KTH - Royal Institute of Technology.

The software comes in two flavors, Python and Matlab scripts to fetch
online data and build the published ecYeast7 GECKO models, and a Python
package which can be used with
`cobrapy <https://opencobra.github.io/cobrapy/>`__ to obtain a ecYeast7
model object, optionally adjusted for provided proteomics data.


Scripts for building a GECKO model
----------------------------------



Required software for the Python module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `Python 2.7 <https://www.python.org/>`__
-  setuptools for python 2.7, accessible
   `here <http://www.lfd.uci.edu/~gohlke/pythonlibs/#setuptools>`__
-  SOAPpy: for this, open command prompt as admin, and then do:

::

    cd C:\Python27\Scripts
    easy_install-2.7 SOAPpy


Required software for the Matlab module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `MATLAB <http://www.mathworks.com/>`__ (7.5 or higher) + Optimization
   Toolbox.
-  The `COBRA toolbox for
   MATLAB <https://github.com/opencobra/cobratoolbox>`__. Note that
   `libSBML <http://sbml.org/Software/libSBML>`__ and the `SBML
   toolbox <http://sbml.org/Software/SBMLToolbox>`__ should both be
   installed. Both of them are free of charge for academic users.
   Aditionally, you should add the cobra folder to your MATLAB search
   path.


Usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See the supporting information of `Sanchez et al.
(2017) <https://dx.doi.org/10.15252/msb.20167411>`__


Using the geckopy Python package for obtaining an adjusted GECKO model object
-----------------------------------------------------------------------------

If all you need is the ecYeast7 model to use together with cobrapy you
can use the ``geckopy`` Python package.

Required software
~~~~~~~~~~~~~~~~~

-  Python 2.7, 3.4, 3.5 or 3.6
-  cobrapy

Installation
~~~~~~~~~~~~

::

    pip install geckopy

.. usage-1:

Usage
~~~~~

.. code:: python

    from geckopy import GeckoModel
    import pandas
    some_measurements = pandas.Series({'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1})
    model = GeckoModel('multi-pool')
    model.limit_proteins(some_measurements)
    model.optimize()
