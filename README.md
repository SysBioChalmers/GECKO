<img src="./GECKO.png" width="200px">

![Current version](https://badge.fury.io/gh/sysbiochalmers%2Fgecko.svg)
[![Gitter](https://badges.gitter.im/SysBioChalmers/GECKO.svg)](https://gitter.im/SysBioChalmers/GECKO)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.7699818.svg)](https://doi.org/10.5281/zenodo.7699818)

## About GECKO 3

The **GECKO** toolbox enhances a **G**enome-scale model to account for **E**nzyme **C**onstraints, using **K**inetics and **O**mics. The resulting enzyme-constrained model (**ecModel**) can be used to perform simulations where enzyme allocation is either drawn from a total protein pool, or constrained by measured protein levels from proteomics data.

💡 In the GECKO folder, `protocol.m` contains instructions on how to reconstruct and analyze an ecModel for _S. cerevisiae_. This demonstrates how many of GECKO's functions can be used.

💡 In the [`GECKO/tutorials`](https://github.com/SysBioChalmers/GECKO/tree/main/tutorials) folder there are examples of how GECKO can be applied to GEMs, in either of its _full_ or _light_ forms. Each `protocol.m` contains instructions on how to reconstruct and analyze an ecModel, demonstrating how different fuctions in GECKO can be used. The source code documentation is also available
[online](http://sysbiochalmers.github.io/GECKO/doc/).

_**Note:** Regarding code and model compatibility with earlier GECKO versions, see [Previous versions: GECKO 1 and 2](#previous-versions-gecko-1-and-2)_.

### Citation

- A GECKO 3 publication is currently under consideration, citation information will appear here in due course.
- For GECKO release 2, please cite [Domenzain et al. (2022) doi:10.1038/s41467-022-31421-1](https://doi.org/10.1038/s41467-022-31421-1).
- For GECKO release 1, please cite [Sánchez et al. (2017) doi:10.15252/msb.20167411](https://doi.org/10.15252/msb.20167411).

### Getting started

#### Required software

- MATLAB version 2019b or later, no additional MathWorks toolboxes are required.
- [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) version 2.7.12 or later. The RAVEN Toolbox Wiki contains [installation instructions for both RAVEN](https://github.com/SysBioChalmers/RAVEN/wiki/Installation) and [Gurobi](https://github.com/SysBioChalmers/RAVEN/wiki/Installation#solvers). Briefly, RAVEN is either downloaded via `git clone`, as ZIP-archive from GitHub, or installed as a [MATLAB AddOn](https://se.mathworks.com/matlabcentral/fileexchange/112330-raven-toolbox). After finishing all installation instructions, the user should run installation checks in MATLAB with: `checkInstallation`.
- [Gurobi Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer/) is recommended for simulations (free academic license available). Alternatively, the open-source [GNU Linear Programming Kit](https://www.gnu.org/software/glpk/) (distributed with RAVEN) or SoPlex as part of the [SCIP Optimization Suite](https://scipopt.org/) can be used.
- [Docker](https://www.docker.com/) for running DLKcat. Installation instructions are available at https://docs.docker.com/get-docker .

#### Installation

- The preferred way to download GECKO is via git clone:

```
git clone --depth=1 https://github.com/SysBioChalmers/GECKO
```

- Alternatively, [a ZIP-archive can be directly downloaded from GitHub](https://github.com/SysBioChalmers/GECKO/releases). The ZIP-archive should be extracted to a disk location where the user has read- and write-access rights.

- After `git clone` or extracting the ZIP-archive, the user should navigate in MATLAB to the GECKO folder. GECKO can then be installed with the command that adds GECKO (sub-)folders to the MATLAB path::

```
  cd('C:\path\to\GECKO') % Modify to match GECKO folder and operating system
  GECKOInstaller.install
```

- If desired, a removal command is available as::

```
  GECKOInstaller.uninstall
```

All set! 🚀

## Previous versions: GECKO 1 and 2

Due to significant refactoring of the code, GECKO version 3 is largely not backwards compatible with earlier GECKO versions.

- Most notably, GECKO 3 ecModels have an `.ec` structure containing all enzyme and kcat information.
- In addition, in GECKO 3 enzymes are incorporated in the S-matrix as MW/kcat, while in GECKO 1 and 2 this was 1/kcat (where the MW was instead considered in the protein exchange reactions).
- GECKO 3 ecModels can be stored in YAML file format that retains all model content.
- Most functions in GECKO 3 do not work on ecModels generated with GECKO versions 1 or 2.
- ecModels generated in GECKO 3 do not work with functions from GECKO versions 1 or 2.
- At this moment, there are no Python functions to work with GECKO 3 formatted ecModels.
- The last GECKO 2 release (2.0.3) is available [here](https://github.com/SysBioChalmers/GECKO/releases/tag/v2.0.3).
- The `gecko2` branch remains available for any potential fixes.

## Contributing

Contributions are always welcome! Please read the [contributing guidelines](https://github.com/SysBioChalmers/GECKO/blob/main/.github/CONTRIBUTING.md) to get started.
