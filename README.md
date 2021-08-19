
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- [![DOI](https://zenodo.org/badge/331070175.svg)](https://zenodo.org/badge/latestdoi/331070175) -->

# An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example

This GitHub repository provides the code of the tutorial on how to
implement time-dependent cohort state-transition models (cSTMs) in R
using a cost-effectiveness analysis (CEA) example, explained in the
following manuscript:

-   Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM,
    Pechlivanoglou P, Jalal H. [An Introductory Tutorial on Cohort
    State-Transition Models in R Using a Cost-Effectiveness Analysis
    Example](http://arxiv.org/abs/2001.07824). arXiv:200107824v3.
    2021:1-26.

The
[`R`](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/tree/main/R)
folder includes two different scripts corresponding to functions used to
synthesize cSTMs outputs and conduct several sensitivity analyses: -
[`Funtions.R`](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/R/Functions.R):
Functions to generate epidemiological measures from time-dependent
cSTMs. -
[`Functions_cSTM_time_indep.R`](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/blob/main/R/Functions_cSTM_time_indep.R):
These functions wrap the time-dependent cSTM, compute CEA measures, and
generate probabilistic sensitivity analysis (PSA) input datasets.

We recommend familiarizing with the [DARTH](http://darthworkgroup.com)
coding framework described in

-   Alarid-Escudero F, Krijkamp EM, Pechlivanoglou P, Jalal HJ, Kao SYZ,
    Yang A, Enns EA. [A Need for Change! A Coding Framework for
    Improving Transparency in Decision
    Modeling](https://link.springer.com/article/10.1007/s40273-019-00837-x).
    [PharmacoEconomics](https://www.springer.com/journal/40273),
    2190;37(11):1329–1339. <https://doi.org/10.1007/s40273-019-00837-x>

To run the CEA, you require [`dampack`: Decision-Analytic Modeling
Package](https://cran.r-project.org/web/packages/dampack/index.html), an
R package for analyzing and visualizing the health economic outputs of
decision models.

## Full list of Contributors:

-   [Fernando Alarid-Escudero](https://github.com/feralaes)

-   [Eline Krijkamp](https://github.com/krijkamp)

-   [Eva Enns](https://github.com/evaenns)

-   [Alan Yang](https://github.com/alanyang0924)

-   [Myriam
    Hunink](http://www.erasmus-epidemiology.nl/people/profile.php?id=45)

-   [Petros Pechlivanoglou](https://github.com/ppehli)

-   [Hawre Jalal](https://github.com/hjalal)
