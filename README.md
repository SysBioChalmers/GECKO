![](GECKO.png?raw=true)
========================

The **GECKO** toolbox is a Matlab/Python package for enhancing a **G**enome-scale model to account for **E**nzyme **C**onstraints, using **K**inetics and **O**mics. It is the companion software to the following publication:

_Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee, Eduard J. Kerkhoven, Jens Nielsen (2017). Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic constraints. Molecular Systems Biology, 13(8): 935_ (http://www.dx.doi.org/10.15252/msb.20167411)

GECKO was programmed by Benjamin J. Sanchez (@BenjaSanchez), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology, and Cheng Zhang, Science for Life Laboratory, KTH - Royal Institute of Technology.

Last update: 2017-09-11

========================

## Required Software for Python module:

* Python 2.7.X (_https://www.python.org/_)
* setuptools for python 2.7.X (_http://www.lfd.uci.edu/~gohlke/pythonlibs/#setuptools_)
* SOAPpy: for this, open command promt as admin, and then type:
    _> cd C:\Python27\Scripts_
    _> easy_install-2.7 SOAPpy_

========================

## Required Software for Matlab module:

* MATLAB (7.5 or higher) + Optimization Toolbox (http://www.mathworks.com/).
* The COBRA toolbox for MATLAB (https://github.com/opencobra/cobratoolbox). Note that libSBML (http://sbml.org/Software/libSBML) and the SBML toolbox (http://sbml.org/Software/SBMLToolbox) should both be installed. Both of them are free of charge for academic users. Aditionally, you should add the cobra folder to your MATLAB search path.

========================