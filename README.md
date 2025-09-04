# ICPM
Integrated population and migratory connectivity model.


[![License:
GPL-2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)



## TODO

- [ ] NIMBLE doesn't like nodes with dimensions of length 1 : how to deal with it (e.g. only 1 cluster in a season) ?


## Content

This project is structured as follow:

```
.
├─ README.md                                  # Presentation of the project
├─ DESCRIPTION                                # Project metadata
├─ LICENSE.md                                 # License of the project
|
├─ pkg_assemble_submodels/                    # Contains a small package to assemble submodels units into a single NIMBLE code
|  ├─ R/                                      # Contains R functions of the package
|  ├─ README                                  # Presentation of the package
|  ├─ make.R                                  # Small script to show how to use the package
|  └─ _targets.R                              # Targets file of the package in order to visualize it
|
├─ R_other/                                   # Contains other R functions (not tidied)
|  └─ ...
|
├─ example/                                   # A project example
|  ├─ workflow.R                              # Workflow
|  └─ submodels/                              # Contains NIMBLE codes of submodels used in this example
|     └─ ...
|
├─ submodels/                                 # Contains examples of ICPM submodels
|  └─ ...
|
├─ simulation/                                # Contains scripts to simulate data and simulated datasets
|  └─ ...
|
├─ test/                                      # R scripts to test some parts of the code
|  └─ ...
|
├─ cmr_pochard/                               # CMR model for common pochard
|  └─ ...
|
├─ .gitignore                                 # files/folders not to track with git 
└─ ICPM.Rproj                                 # Rstudio project
```

## Package to assemble submodels

The directory *pkg_assemble_sumbodels/* contains an independant R package. The aim of this very small package is to automatically create a NIMBLE code for an integrated model, given separate codes for submodels and priors.
The package concatenates codes from selected submodels, scans for parameters that need priors and add priors provided by the user.

To load this package, use :
`devtools::load_all(<path_to_pkg>)`
for example `devtools::load_all('pkg_assemble_submodels')` if you are working in the root directory of this compendium.

## Data

The package can be used to run simulations, in which case external data is not needed.
The model can also be applied to external data.
> [!NOTE] TODO
> Documentation on the formats needed for the data.


## Installation

To install this compendium:

- [Fork](https://docs.github.com/en/get-started/quickstart/contributing-to-projects)
  this [repository](https://github.com/cdoucerain/ICPM_simulations.git) using the GitHub interface.
- Open [RStudio IDE](https://posit.co/products/open-source/rstudio/) and create a 
  **New Project** from **Version Control** to [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
  your fork.



## License

The project is released under the GPL-2 license. 
See LICENSE.md for description of the license (https://choosealicense.com/licenses/gpl-2.0/)

