
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![check](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/actions/workflows/check.yaml/badge.svg)](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/actions/workflows/check.yaml)
[![DOI](https://zenodo.org/badge/357362984.svg)](https://zenodo.org/badge/latestdoi/357362984)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

# An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example

This GitHub repository has the R code for a tutorial on implementing
time-independent cohort state-transition models (cSTMs) in R, using a
cost-effectiveness analysis (CEA) example. The tutorial is the following
article:

- Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
  Pechlivanoglou P, Jalal H. [An Introductory Tutorial on Cohort
  State-Transition Models in R Using a Cost-Effectiveness Analysis
  Example](https://journals.sagepub.com/doi/full/10.1177/0272989X221103163).
  [Medical Decision Making](https://journals.sagepub.com/home/mdm),
  2023;43(1):3-20. <https://doi.org/10.1177/0272989X221103163>

Every release of the code is archived in Zenodo:
<https://doi.org/10.5281/zenodo.5223092> (this DOI always resolves to
the latest version).

> **Correction (v0.3.0):** Versions up to v0.2.1, and the published
> article, assigned Simpson’s 1/3rd within-cycle correction weights by
> vector position instead of by time. They used `1/3, 2/3, 4/3, …`
> instead of `1/3, 4/3, 2/3, …`. v0.3.0 fixes this. With the fix,
> base-case QALYs are about 0.32 higher per strategy and ICERs change by
> about 0.4%; the conclusions are unchanged. See [`NEWS.md`](NEWS.md)
> for details and a before/after table. We thank Lauren N. Steimle and
> [ambtnhlr](https://github.com/ambtnhlr) for finding and reporting this
> error (see [Acknowledgements](#acknowledgements)).

## Quick start

``` r
# 1. Clone the repository (or download it as a ZIP from GitHub) and open
#    `cohort-modeling-tutorial-intro.Rproj` in RStudio, or set the working
#    directory to the repository's root folder.

# 2. Install the exact package versions used by the tutorial
install.packages("renv")
renv::restore()

# 3. Run the tutorial: base case, CEA, PSA and value of information
source("analysis/cSTM_time_indep.R")

# 4. (Optional) Run the unit tests
testthat::test_dir("tests/testthat")
```

The tutorial was most recently checked with R 4.5.3. Package versions
are pinned in [`renv.lock`](renv.lock). If you would rather not use
`renv`, install the packages listed under `Imports` in
[`DESCRIPTION`](DESCRIPTION) with `install.packages()`.

## Repository structure

    ├── analysis/
    │   └── cSTM_time_indep.R            # Main script: Sick-Sicker cSTM, CEA, PSA, VOI
    ├── R/
    │   ├── Functions.R                  # Epidemiological outputs, within-cycle correction,
    │   │                                # rate-to-probability conversion, matrix checks,
    │   │                                # and CEA/PSA functions (copies from dampack)
    │   └── Functions_cSTM_time_indep.R  # Model wrapper, CEA outputs and PSA input generation
    ├── tests/testthat/                  # Unit and regression tests
    ├── manuscript/                      # Archival sources of the published article (see manuscript/README.md)
    ├── DESCRIPTION                      # Project version and dependencies
    ├── renv.lock                        # Pinned package versions
    ├── NEWS.md                          # Changelog
    └── CITATION.cff                     # Citation metadata

The `R` folder includes two scripts of functions used to synthesize cSTM
outputs and run several sensitivity analyses:

- [`Functions.R`](R/Functions.R): Functions that generate
  epidemiological measures from time-independent cSTMs and compute the
  within-cycle correction, parameter transformations, matrix checks, and
  CEA and PSA visualizations.
- [`Functions_cSTM_time_indep.R`](R/Functions_cSTM_time_indep.R):
  Functions that wrap the time-independent cSTM, compute CEA measures,
  and generate probabilistic sensitivity analysis (PSA) input datasets.

## Dependencies on `dampack` and `darthtools`

The repository is self-contained: `R/Functions.R` includes local copies
of the functions it uses from two DARTH packages:

- [`dampack`](https://cran.r-project.org/package=dampack):
  Decision-Analytic Modeling Package, an R package for analyzing and
  visualizing the health economic outputs of decision models.
- [`darthtools`](https://github.com/DARTH-git/darthtools): helper
  functions for the within-cycle correction, rate conversion and
  transition-matrix checks.

You don’t need to install either package to run the tutorial. To use
them anyway:

``` r
install.packages("dampack")                          # from CRAN
# install.packages("remotes")
remotes::install_github("DARTH-git/darthtools")      # from GitHub
```

**Note:** `darthtools` up to v0.4.0 has the same Simpson’s 1/3rd
indexing error as this repository had. Use the corrected `gen_wcc()` in
`R/Functions.R` until a fixed `darthtools` is released.

We recommend becoming familiar with the
[DARTH](https://darthworkgroup.com) coding framework described in

- Alarid-Escudero F, Krijkamp EM, Pechlivanoglou P, Jalal HJ, Kao SYZ,
  Yang A, Enns EA. [A Need for Change! A Coding Framework for Improving
  Transparency in Decision
  Modeling](https://link.springer.com/article/10.1007/s40273-019-00837-x).
  [PharmacoEconomics](https://www.springer.com/journal/40273),
  2019;37(11):1329–1339. <https://doi.org/10.1007/s40273-019-00837-x>

## How to cite this R code in your article

You can cite the R code in this repository like this: “we based our
analysis on the R code from Alarid-Escudero F et al. (2023)”. Include
both of the following references in your reference list. Update the
version, version DOI and ‘last accessed’ date to match the release you
used. Each version’s DOI is listed on the [Zenodo
record](https://doi.org/10.5281/zenodo.5223092), and the GitHub “Cite
this repository” button gives the metadata from
[`CITATION.cff`](CITATION.cff).

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
> Pechlivanoglou P, Jalal H. An Introductory Tutorial on Cohort
> State-Transition Models in R Using a Cost-Effectiveness Analysis
> Example. Medical Decision Making, 2023;43(1):3-20.
> <https://doi.org/10.1177/0272989X221103163>

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
> Pechlivanoglou P, Jalal H. R Code for An Introductory Tutorial on
> Cohort State-Transition Models in R Using a Cost-Effectiveness
> Analysis Example (Version v0.3.0). Zenodo.
> [10.5281/zenodo.5223092](https://doi.org/10.5281/zenodo.5223092). Last
> accessed DD Month YYYY.

If you adapted the code, write “Adapted from:” or “Based on” so readers
know you modified it. For more on how and when to cite computer code,
see [Writing Code (from the MIT Academic Integrity
Handbook)](https://integrity.mit.edu/handbook/writing-code).

## Reporting issues and contributing

Found a bug or have a suggestion? Please [open an
issue](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues).
Include a minimal reproducible example and the output of
`sessionInfo()`.

Pull requests are welcome. Before you submit one, please:

1.  Add or update tests in `tests/testthat/` that cover the change.
2.  Make sure `testthat::test_dir("tests/testthat")` passes and that
    `analysis/cSTM_time_indep.R` runs end to end. Continuous integration
    checks both on every pull request.
3.  If you add or update packages, run `renv::snapshot()`.
4.  Describe the change in [`NEWS.md`](NEWS.md).

## Advanced state-transition cohort modeling

To learn how to implement time-dependent cohort state-transition models
(cSTMs) in R using a cost-effectiveness analysis (CEA) example, see the
following article:

- Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
  Pechlivanoglou P, Jalal H. [A Tutorial on Time-Dependent Cohort
  State-Transition Models in R using a Cost-Effectiveness Analysis
  Example](https://journals.sagepub.com/doi/full/10.1177/0272989X221121747).
  [Medical Decision Making](https://journals.sagepub.com/home/mdm).
  2023;43(1):21-41. <https://doi.org/10.1177/0272989X221121747>

## Acknowledgements

We are grateful to the users who have helped improve this code:

- **Lauren N. Steimle** (Georgia Institute of Technology) and
  **[ambtnhlr](https://github.com/ambtnhlr)** each independently found
  and reported the indexing error in the Simpson’s 1/3rd within-cycle
  correction (`gen_wcc()`) and proposed a correction. Their reports led
  to the fix in v0.3.0
  ([\#2](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues/2)).

If you find an error, please [let us
know](#reporting-issues-and-contributing).

## Full list of Contributors:

- [Fernando Alarid-Escudero](https://github.com/feralaes)

- [Eline Krijkamp](https://github.com/krijkamp)

- [Eva Enns](https://github.com/evaenns)

- [Alan Yang](https://github.com/alanyang0924)

- [Myriam
  Hunink](http://www.erasmus-epidemiology.nl/people/profile.php?id=45)

- [Petros Pechlivanoglou](https://github.com/ppehli)

- [Hawre Jalal](https://github.com/hjalal)
