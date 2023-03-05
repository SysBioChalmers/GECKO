.. image:: GECKO.png
   :align: center

|Current Version| |Tests passing| |Gitter| |Zenodo|

About GECKO
-----------

The **GECKO** toolbox is able to enhance a **G**\ enome-scale model to account
for **E**\ nzyme **C**\ onstraints, using **K**\ inetics and **O**\ mics. The resulting enzyme-constrained model
(**ecModel**) can be used to perform simulations where enzyme allocation is either drawn from a total protein pool, or
constrained by measured protein levels from proteomics data.

**Note:** Due to significant refactoring of the code, ecModels generated with GECKO versions 1 or 2 are not compatible
with GECKO 3, and *vice versa*. The latest GECKO 2 release is available `here <https://github.com/SysBioChalmers/GECKO/releases/tag/v2.0.3>`_,
while the ``gecko2`` branch is retained.

**Citation**

- A GECKO 3 publication is currently under consideration, citation information will appear here in due course.
- For GECKO release 2, please cite `Domenzain et al. (2022) <https://doi.org/10.1038/s41467-022-31421-1>`_.
- For GECKO release 1, please cite `Sánchez et al. (2017) <https://doi.org/10.15252/msb.20167411>`_.

Last update: 2023-03-05

Required software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-	MATLAB version 2019b or later, no additional MathWorks toolboxes are required.
-	`RAVEN <https://github.com/SysBioChalmers/RAVEN>`_ Toolbox version 2.7.12 or later.
-	`Gurobi Optimizer <https://www.gurobi.com/solutions/gurobi-optimizer/>`_ is recommended for simulations (free academic license available). Alternatively, the open-source GNU Linear Programming Kit (`GLPK <https://www.gnu.org/software/glpk/>`_, distributed with RAVEN) or SoPlex as part of the `SCIP Optimization Suite <https://scipopt.org/>`_ can be used.
-	`Docker <https://www.docker.com/>`_ for running DLKcat.

Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**GECKO toolbox**

- The preferred way to download GECKO is via git clone::

   git clone --depth=1 https://github.com/SysBioChalmers/GECKO

- Alternatively, a `ZIP-archive <https://github.com/SysBioChalmers/GECKO/releases>`_ can be directly downloaded from GitHub. The ZIP-archive should be extracted to a disk location where the user has read- and write-access rights.

- After git clone or extracting the ZIP-archive, the user should navigate in MATLAB to the GECKO folder. GECKO can then be installed with the command that adds GECKO (sub-)folders to the MATLAB path::

   cd('C:\path\to\GECKO') % Modify to match GECKO folder and OS
   GECKOInstaller.install

- If desired, a removal command is available as::

   GECKOInstaller.uninstall

**RAVEN Toolbox and Gurobi**

- The RAVEN Toolbox Wiki contains installation instructions for both `RAVEN Toolbox <https://github.com/SysBioChalmers/RAVEN/wiki/Installation>`_ and `Gurobi <https://github.com/SysBioChalmers/RAVEN/wiki/Installation#solvers>`_. 

  - Briefly, RAVEN is either downloaded via git clone, as ZIP-archive from GitHub, or installed as `MATLAB AddOn <https://se.mathworks.com/matlabcentral/fileexchange/112330-raven-toolbox>`_.

- After finishing all installation instructions, the user should run installation checks in MATLAB with::

   checkInstallation

**Docker**

- Installation instructions are available at https://docs.docker.com/get-docker/.

Getting started
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the GECKO folder, ``protocols.m`` contains instructions on how to reconstruct 
and analyze an ecModel for *S. cerevisiae*.

Contributing
------------

Contributions are always welcome! Please read the `contributing guidelines <https://github.com/SysBioChalmers/GECKO/blob/devel/.github/CONTRIBUTING.md>`_ to get started.

.. |Current Version| image:: https://badge.fury.io/gh/sysbiochalmers%2Fgecko.svg
   :target: https://badge.fury.io/gh/sysbiochalmers%2Fgecko
.. |Tests passing| image:: https://github.com/SysBioChalmers/GECKO/actions/workflows/tests.yml/badge.svg?branch=main
   :target: https://github.com/SysBioChalmers/GECKO/actions
.. |Gitter| image:: https://badges.gitter.im/SysBioChalmers/GECKO.svg
   :alt: Join the chat at https://gitter.im/SysBioChalmers/GECKO
   :target: https://gitter.im/SysBioChalmers/GECKO?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.7699818.svg
   :target: https://doi.org/10.5281/zenodo.7699818
