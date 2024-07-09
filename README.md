<img src="./GECKO.png" width="200px">

[![Current release](https://img.shields.io/github/release/SysBioChalmers/GECKO/all.svg)](https://GitHub.com/SysBioChalmers/GECKO/releases/)
[![GitHub Discussions](https://img.shields.io/github/discussions-search?query=repo%3Asysbiochalmers%2Fgecko&label=GitHub%20Discussions)](https://github.com/SysBioChalmers/GECKO/discussions)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.7699818.svg)](https://doi.org/10.5281/zenodo.7699818)
[![MATLAB File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/125960-gecko-toolbox)


## About GECKO 3

The **GECKO** toolbox enhances a **G**enome-scale model to account for **E**nzyme **C**onstraints, using **K**inetics and **O**mics. The resulting enzyme-constrained model (**ecModel**) can be used to perform simulations where enzyme allocation is either drawn from a total protein pool, or constrained by measured protein levels from proteomics data.

💡 In the [`GECKO/tutorials`](https://github.com/SysBioChalmers/GECKO/tree/main/tutorials) folder there are examples of how GECKO can be applied to GEMs, in either of its _full_ or _light_ forms. Each `protocol.m` contains instructions on how to reconstruct and analyze an ecModel, demonstrating how different fuctions in GECKO can be used. These two scripts complement the [Nature Protocols](https://doi.org/10.1038/s41596-023-00931-7) paper ([**PDF**](https://drive.google.com/file/d/1_AGz6GmQPOyshfUZ6K-L2myzU43rQ6tu/view)).

### Significant changes since protocol publication
- GECKO **3.2.0**: all protein usage reactions draw from the protein pool, even if they are constrained by proteomics data. This affects **Step 58** in the protocol, changing behaviour of `constrainEnzConcs` and making `updateProtPool` obsolete, `tutorials/full_ecModel/protocol.m` is updated to reflect this change. See [#357](https://github.com/SysBioChalmers/GECKO/issues/375) for more details.  
  
_**Note:** Regarding code and model compatibility with earlier GECKO versions, see [Previous versions: GECKO 1 and 2](https://github.com/SysBioChalmers/GECKO/wiki/Previous-versions:-GECKO-1-and-2)_.

## Cite us

If you have used GECKO in your work, please cite:

> Chen Y, Gustafsson J, Tafur Rangel A, Anton M, Domenzain I, Kittikunapong C, Li F, Yuan L, Nielsen J & Kerkhoven EJ. Reconstruction, simulation and analysis of enzyme-constrained metabolic models using GECKO Toolbox 3.0. *Nature Protocols* **19**, 629-667 (2024). doi:[10.1038/s41596-023-00931-7](https://doi.org/10.1038/s41596-023-00931-7).

## Documentation
**Installation instructions** and other information is available in the **[Wiki](https://github.com/SysBioChalmers/GECKO/wiki)**. The source code documentation is also available 
[online](http://sysbiochalmers.github.io/GECKO/doc/). Use [GitHub Discussions](https://github.com/SysBioChalmers/GECKO/discussions) for support, to ask questions or leave comments.

## Contributing

Contributions are always welcome! Please read the [contributing guidelines](https://github.com/SysBioChalmers/GECKO/blob/main/.github/CONTRIBUTING.md) to get started.
