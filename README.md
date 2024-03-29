
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/357362984.svg)](https://zenodo.org/badge/latestdoi/357362984)

# An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example

This GitHub repository provides the code of the tutorial on how to
implement time-independent cohort state-transition models (cSTMs) in R
using a cost-effectiveness analysis (CEA) example, explained in the
following manuscript:

- Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
  Pechlivanoglou P, Jalal H. [An Introductory Tutorial on Cohort
  State-Transition Models in R Using a Cost-Effectiveness Analysis
  Example](https://journals.sagepub.com/doi/full/10.1177/0272989X221103163).
  [Medical Decision Making](https://journals.sagepub.com/home/mdm),
  2023;43(1):3-20. <https://doi.org/10.1177/0272989X221103163>

The release that accompanies the published article has been archived in
zenodo: <https://zenodo.org/badge/latestdoi/357362984>

The
[`R`](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/tree/main/R)
folder includes two different scripts corresponding to functions used to
synthesize cSTMs outputs and conduct several sensitivity analyses:

- [`Funtions.R`](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/R/Functions.R):
  Functions that generate epidemiological measures from time-independent
  cSTMs and compute within-cycle correction, parameter transformation,
  matrix checks, and CEA and PSA visualization.
- [`Functions_cSTM_time_indep.R`](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/R/Functions_cSTM_time_indep.R):
  These functions wrap the time-independent cSTM, compute CEA measures,
  and generate probabilistic sensitivity analysis (PSA) input datasets.

## How to cite this R code in your article

You can cite the R code in this repository like this “we based our
analysis using the R code from Alarid-Escudero F et al. (2022)”. Here is
the full bibliographic reference to include in your reference list for
the manuscript and the R code (don’t forget to update the ‘last
accessed’ date):

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
> Pechlivanoglou P, Jalal H. An Introductory Tutorial on Cohort
> State-Transition Models in R Using a Cost-Effectiveness Analysis
> Example. Medical Decision Making, 2023;43(1):3-20.

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
> Pechlivanoglou P, Jalal H (2022). R Code for An Introductory Tutorial
> on Cohort State-Transition Models in R Using a Cost-Effectiveness
> Analysis Example (Version v0.2.1). Zenodo.
> [10.5281/zenodo.5223093](https://www.doi.org/10.5281/zenodo.5223093).
> Last accessed 30 March 2022.

If you adapted the code, you should indicate “Adapted from:” or “Based
on” so it is understood that you modified the code. For more information
on how to cite computer code, we refer the user to review [Writing Code
(from MIT Research
Guide)](https://integrity.mit.edu/handbook/writing-code), which provides
examples of how and when to cite computer code.

## Advanced state-transition cohort modeling

To learn more on on how to implement time-dependent cohort
state-transition models (cSTMs) in R using a cost-effectiveness analysis
(CEA) example, we refer the reader to the following manuscript:

- Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
  Pechlivanoglou P, Jalal H. [A Tutorial on Time-Dependent Cohort
  State-Transition Models in R using a Cost-Effectiveness Analysis
  Example](https://journals.sagepub.com/doi/full/10.1177/0272989X221121747).
  [Medical Decision Making](https://journals.sagepub.com/home/mdm).
  2023;43(1):21-41. <https://doi.org/10.1177/0272989X221121747>

## Preliminaries

- Install [RStudio](https://www.rstudio.com/products/rstudio/download/)
- Install
  [`dampack`](https://cran.r-project.org/web/packages/dampack/index.html)
  R package from CRAN

``` r
# Install release version from CRAN
install.packages("dampack")

# Or install development version from GitHub
# devtools::install_github("DARTH-git/dampack")
```

- Install `devtools` to install
  [`darthtools`](https://github.com/DARTH-git/darthtools) R package from
  [DARTH’s GitHub](https://github.com/DARTH-git)

``` r
# Install release version from CRAN
install.packages("devtools")

# Or install development version from GitHub
# devtools::install_github("r-lib/devtools")
```

- Install `darthtools` using `devtools`

``` r
# Install development version from GitHub
devtools::install_github("DARTH-git/darthtools")
```

We recommend familiarizing with the [DARTH](http://darthworkgroup.com)
coding framework described in

- Alarid-Escudero F, Krijkamp EM, Pechlivanoglou P, Jalal HJ, Kao SYZ,
  Yang A, Enns EA. [A Need for Change! A Coding Framework for Improving
  Transparency in Decision
  Modeling](https://link.springer.com/article/10.1007/s40273-019-00837-x).
  [PharmacoEconomics](https://www.springer.com/journal/40273),
  2190;37(11):1329–1339. <https://doi.org/10.1007/s40273-019-00837-x>

To run the CEA, you require [`dampack`: Decision-Analytic Modeling
Package](https://cran.r-project.org/web/packages/dampack/index.html), an
R package for analyzing and visualizing the health economic outputs of
decision models.

## Use repository as a regular coding template

1.  On the [tutorial’s GitHub
    repository](https://github.com/DARTH-git/cohort-modeling-tutorial-intro),
    navigate to the main page of the repository
    (<https://github.com/DARTH-git/cohort-modeling-tutorial-intro>).
2.  Above the file list, click **Clone or download** and select either
    1.  **Open in desktop**, which requires the user to have a GitHub
        desktop installed, or
    2.  **Download zip** that will ask the user to download the whole
        repository as a .zip file.
3.  Open the RStudio project `cohort-modeling-tutorial-intro.Rproj`.
4.  Install all the required packages (as mentioned above)
    - [`dampack`](https://cran.r-project.org/web/packages/dampack/index.html)
    - [`darthtools`](https://github.com/DARTH-git/darthtools)
5.  Run the scripts in the analysis folder.
6.  Modify or adapt these scripts as needed for your project or
    analysis.

## Full list of Contributors:

- [Fernando Alarid-Escudero](https://github.com/feralaes)

- [Eline Krijkamp](https://github.com/krijkamp)

- [Eva Enns](https://github.com/evaenns)

- [Alan Yang](https://github.com/alanyang0924)

- [Myriam
  Hunink](http://www.erasmus-epidemiology.nl/people/profile.php?id=45)

- [Petros Pechlivanoglou](https://github.com/ppehli)

- [Hawre Jalal](https://github.com/hjalal)
